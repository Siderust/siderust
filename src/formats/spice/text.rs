// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Generic parser for SPICE text kernels.
//!
//! ## Scientific scope
//!
//! This module is format-oriented. It parses text-kernel assignments for frame,
//! leapsecond, instrument, and body-orientation metadata without interpreting
//! the scientific meaning of most keys.
//!
//! ## Technical scope
//!
//! Supported syntax covers `\begindata` / `\begintext` blocks, scalar values,
//! parenthesised arrays, quoted strings, and `+=` array-append assignments.
//!
//! ## References
//!
//! - NAIF. *Kernel Required Reading*.
//! - NAIF. *Text Kernel Required Reading*.

use std::collections::HashMap;

use super::SpiceError;

/// Value parsed from a SPICE text kernel.
#[derive(Debug, Clone, PartialEq)]
pub enum TextValue {
    /// Signed integer scalar.
    Integer(i64),
    /// Floating-point scalar.
    Float(f64),
    /// Text scalar, either quoted or a non-numeric bare token.
    Text(String),
    /// Parenthesised array of values.
    Array(Vec<TextValue>),
}

/// Parsed SPICE text kernel: a map from uppercase key to value.
///
/// # Examples
///
/// ```rust
/// use siderust::formats::spice::{TextKernel, text::TextValue};
///
/// let kernel = TextKernel::parse("\\begindata\nANSWER = 42\n")?;
/// assert_eq!(kernel.get("answer"), Some(&TextValue::Integer(42)));
/// # Ok::<_, siderust::formats::spice::SpiceError>(())
/// ```
#[derive(Debug, Default, Clone)]
pub struct TextKernel {
    /// Parsed uppercase assignment keys and values.
    pub data: HashMap<String, TextValue>,
}

impl TextKernel {
    /// Parse a text kernel from a string slice.
    pub fn parse(src: &str) -> Result<Self, SpiceError> {
        let mut kernel = Self::default();
        let mut in_data = false;
        let mut statement = String::new();

        for raw_line in src.lines() {
            let trimmed = raw_line.trim();
            if trimmed.eq_ignore_ascii_case("\\begindata") {
                in_data = true;
                statement.clear();
                continue;
            }
            if trimmed.eq_ignore_ascii_case("\\begintext") {
                in_data = false;
                statement.clear();
                continue;
            }
            if !in_data {
                continue;
            }
            if trimmed.starts_with("/*") && trimmed.ends_with("*/") {
                continue;
            }

            let clean = strip_comment_fragments(raw_line);
            let clean = clean.trim();
            if clean.is_empty() {
                continue;
            }

            if !statement.is_empty() {
                statement.push(' ');
            }
            statement.push_str(clean);

            if !statement_complete(&statement) {
                continue;
            }

            let (key, append, value) = parse_assignment(&statement)?;
            store_value(&mut kernel.data, key, append, value);
            statement.clear();
        }

        if !statement.trim().is_empty() {
            return Err(SpiceError::FormatParse(
                "unterminated text-kernel assignment".to_string(),
            ));
        }

        Ok(kernel)
    }

    /// Retrieve a key's value, or `None` if absent.
    pub fn get(&self, key: &str) -> Option<&TextValue> {
        self.data.get(&key.to_ascii_uppercase())
    }

    /// Retrieve a key as an `f64` array, flattening nested arrays one level.
    pub fn get_f64_array(&self, key: &str) -> Option<Vec<f64>> {
        match self.get(key)? {
            TextValue::Integer(value) => Some(vec![*value as f64]),
            TextValue::Float(value) => Some(vec![*value]),
            TextValue::Array(values) => flatten_f64(values),
            TextValue::Text(_) => None,
        }
    }

    /// Retrieve a key as an `i64` array, flattening nested arrays one level.
    pub fn get_i64_array(&self, key: &str) -> Option<Vec<i64>> {
        match self.get(key)? {
            TextValue::Integer(value) => Some(vec![*value]),
            TextValue::Array(values) => flatten_i64(values),
            TextValue::Float(_) | TextValue::Text(_) => None,
        }
    }
}

fn strip_comment_fragments(line: &str) -> String {
    let mut out = String::new();
    let mut rest = line;
    while let Some(start) = rest.find("/*") {
        out.push_str(&rest[..start]);
        let Some(end) = rest[start + 2..].find("*/") else {
            return out;
        };
        rest = &rest[start + 2 + end + 2..];
    }
    out.push_str(rest);
    out
}

fn statement_complete(statement: &str) -> bool {
    let mut paren_depth = 0_i32;
    let mut quote = None;
    for ch in statement.chars() {
        if let Some(active) = quote {
            if ch == active {
                quote = None;
            }
            continue;
        }
        match ch {
            '\'' | '"' => quote = Some(ch),
            '(' => paren_depth += 1,
            ')' => paren_depth -= 1,
            _ => {}
        }
    }
    quote.is_none() && paren_depth == 0
}

fn parse_assignment(statement: &str) -> Result<(String, bool, TextValue), SpiceError> {
    let (key, value, append) = if let Some(idx) = statement.find("+=") {
        (&statement[..idx], &statement[idx + 2..], true)
    } else if let Some(idx) = statement.find('=') {
        (&statement[..idx], &statement[idx + 1..], false)
    } else {
        return Err(SpiceError::FormatParse(format!(
            "text kernel assignment missing '=': {statement}"
        )));
    };

    let key = key.trim();
    if key.is_empty() {
        return Err(SpiceError::FormatParse(format!(
            "text kernel assignment missing key: {statement}"
        )));
    }

    let mut parser = ValueParser::new(value.trim());
    let parsed = parser.parse_value()?;
    parser.skip_delimiters();
    if !parser.is_eof() {
        return Err(SpiceError::FormatParse(format!(
            "unexpected trailing text in assignment '{statement}'"
        )));
    }

    Ok((key.to_ascii_uppercase(), append, parsed))
}

fn store_value(data: &mut HashMap<String, TextValue>, key: String, append: bool, value: TextValue) {
    if !append {
        data.insert(key, value);
        return;
    }

    match data.remove(&key) {
        Some(TextValue::Array(mut existing)) => {
            match value {
                TextValue::Array(mut extra) => existing.append(&mut extra),
                other => existing.push(other),
            }
            data.insert(key, TextValue::Array(existing));
        }
        Some(existing) => {
            let mut values = vec![existing];
            match value {
                TextValue::Array(mut extra) => values.append(&mut extra),
                other => values.push(other),
            }
            data.insert(key, TextValue::Array(values));
        }
        None => {
            let value = match value {
                TextValue::Array(values) => TextValue::Array(values),
                other => TextValue::Array(vec![other]),
            };
            data.insert(key, value);
        }
    }
}

fn flatten_f64(values: &[TextValue]) -> Option<Vec<f64>> {
    let mut out = Vec::new();
    for value in values {
        match value {
            TextValue::Integer(number) => out.push(*number as f64),
            TextValue::Float(number) => out.push(*number),
            TextValue::Array(inner) => {
                for nested in inner {
                    match nested {
                        TextValue::Integer(number) => out.push(*number as f64),
                        TextValue::Float(number) => out.push(*number),
                        TextValue::Text(_) | TextValue::Array(_) => return None,
                    }
                }
            }
            TextValue::Text(_) => return None,
        }
    }
    Some(out)
}

fn flatten_i64(values: &[TextValue]) -> Option<Vec<i64>> {
    let mut out = Vec::new();
    for value in values {
        match value {
            TextValue::Integer(number) => out.push(*number),
            TextValue::Array(inner) => {
                for nested in inner {
                    match nested {
                        TextValue::Integer(number) => out.push(*number),
                        TextValue::Float(_) | TextValue::Text(_) | TextValue::Array(_) => {
                            return None
                        }
                    }
                }
            }
            TextValue::Float(_) | TextValue::Text(_) => return None,
        }
    }
    Some(out)
}

struct ValueParser<'a> {
    src: &'a str,
    index: usize,
}

impl<'a> ValueParser<'a> {
    fn new(src: &'a str) -> Self {
        Self { src, index: 0 }
    }

    fn parse_value(&mut self) -> Result<TextValue, SpiceError> {
        self.skip_delimiters();
        let Some(ch) = self.peek() else {
            return Err(SpiceError::FormatParse(
                "missing text-kernel value".to_string(),
            ));
        };
        match ch {
            '(' => self.parse_array(),
            '\'' | '"' => self.parse_string(),
            _ => self.parse_bare(),
        }
    }

    fn parse_array(&mut self) -> Result<TextValue, SpiceError> {
        self.bump();
        let mut values = Vec::new();
        loop {
            self.skip_delimiters();
            match self.peek() {
                Some(')') => {
                    self.bump();
                    break;
                }
                Some(_) => values.push(self.parse_value()?),
                None => {
                    return Err(SpiceError::FormatParse(
                        "unterminated text-kernel array".to_string(),
                    ));
                }
            }
        }
        Ok(TextValue::Array(values))
    }

    fn parse_string(&mut self) -> Result<TextValue, SpiceError> {
        let quote = self
            .bump()
            .ok_or_else(|| SpiceError::FormatParse("missing quote delimiter".to_string()))?;
        let mut out = String::new();
        while let Some(ch) = self.bump() {
            if ch == quote {
                return Ok(TextValue::Text(out));
            }
            out.push(ch);
        }
        Err(SpiceError::FormatParse(
            "unterminated quoted string in text kernel".to_string(),
        ))
    }

    fn parse_bare(&mut self) -> Result<TextValue, SpiceError> {
        let start = self.index;
        while let Some(ch) = self.peek() {
            if ch.is_whitespace() || ch == ',' || ch == ')' {
                break;
            }
            self.bump();
        }
        let token = self.src[start..self.index].trim();
        if token.is_empty() {
            return Err(SpiceError::FormatParse(
                "empty text-kernel token".to_string(),
            ));
        }
        Ok(parse_scalar(token))
    }

    fn skip_delimiters(&mut self) {
        while let Some(ch) = self.peek() {
            if ch.is_whitespace() || ch == ',' {
                self.bump();
            } else {
                break;
            }
        }
    }

    fn peek(&self) -> Option<char> {
        self.src[self.index..].chars().next()
    }

    fn bump(&mut self) -> Option<char> {
        let ch = self.peek()?;
        self.index += ch.len_utf8();
        Some(ch)
    }

    fn is_eof(&self) -> bool {
        self.index >= self.src.len()
    }
}

fn parse_scalar(token: &str) -> TextValue {
    if let Ok(value) = token.parse::<i64>() {
        return TextValue::Integer(value);
    }
    if let Ok(value) = token.parse::<f64>() {
        return TextValue::Float(value);
    }
    TextValue::Text(token.to_string())
}

#[cfg(test)]
mod tests {
    use super::{TextKernel, TextValue};

    #[test]
    fn parse_empty_text_kernel() {
        let kernel = TextKernel::parse("ignored text only").unwrap();
        assert!(kernel.data.is_empty());
    }

    #[test]
    fn parse_single_integer_key() {
        let kernel = TextKernel::parse("\\begindata\nANSWER = 42\n").unwrap();
        assert_eq!(kernel.get("answer"), Some(&TextValue::Integer(42)));
    }

    #[test]
    fn parse_float_array() {
        let kernel = TextKernel::parse("\\begindata\nVALUES = ( 1.0 2.5 3.0E+01 )\n").unwrap();
        assert_eq!(kernel.get_f64_array("values"), Some(vec![1.0, 2.5, 30.0]));
    }

    #[test]
    fn parse_quoted_string() {
        let kernel = TextKernel::parse("\\begindata\nNAME = 'IAU_EARTH'\n").unwrap();
        assert_eq!(
            kernel.get("NAME"),
            Some(&TextValue::Text("IAU_EARTH".to_string()))
        );
    }

    #[test]
    fn parse_begin_end_data_blocks_with_leading_text() {
        let src = "header\n\\begindata\nX = 1\n\\begintext\ntrailer\n";
        let kernel = TextKernel::parse(src).unwrap();
        assert_eq!(kernel.get("x"), Some(&TextValue::Integer(1)));
    }

    #[test]
    fn parse_plus_equals_array_append() {
        let src = "\\begindata\nA = ( 1 2 )\nA += ( 3 4 )\n";
        let kernel = TextKernel::parse(src).unwrap();
        assert_eq!(kernel.get_i64_array("a"), Some(vec![1, 2, 3, 4]));
    }
}
