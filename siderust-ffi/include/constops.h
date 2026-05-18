/* SPDX-License-Identifier: AGPL-3.0-or-later
 * Copyright (C) 2026 Vallés Puig, Ramon
 *
 * Hand-written C ABI for the constops planning loop.
 *
 * Stable as of siderust-ffi 0.4.x — gated behind the `constops` Cargo feature.
 * The wire format for every entry point is JSON (UTF-8); the request/response
 * shapes mirror the corresponding axum handlers in `constops::api` exactly.
 */

#ifndef CONSTOPS_FFI_H
#define CONSTOPS_FFI_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Status codes returned by every constops_* entry point. */
enum constops_error_code_t
#ifdef __cplusplus
  : int32_t
#endif
{
  CONSTOPS_OK       = 0,
  CONSTOPS_JSON     = 1,
  CONSTOPS_PASS     = 2,
  CONSTOPS_SCORE    = 3,
  CONSTOPS_SCHEDULE = 4,
  CONSTOPS_VALIDATE = 5,
  CONSTOPS_INGEST   = 6,
  CONSTOPS_INTERNAL = 99
};
#ifndef __cplusplus
typedef int32_t constops_error_code_t;
#endif

/*
 * Heap-allocated UTF-8 byte buffer returned by constops_* entry points.
 *
 * `data` is owned by the FFI allocator until released via
 * constops_buffer_free. A zero-length buffer uses `data == NULL` and
 * `len == 0`.
 */
typedef struct ConstopsBuffer {
  uint8_t* data;
  size_t   len;
} ConstopsBuffer;

/*
 * Generate satellite passes from a `GeneratePassesRequest` JSON payload.
 *
 * Returns 0 on success and writes a JSON-encoded `GeneratePassesResponse`
 * into `*out`. Returns a non-zero status on failure; call
 * `constops_last_error()` for a human-readable diagnostic.
 *
 * Caller must release `*out` via `constops_buffer_free` on success.
 */
int32_t constops_passes_generate(const uint8_t* req, size_t req_len, ConstopsBuffer* out);

/* Score link opportunities from a `ScoreLinksRequest` JSON payload. */
int32_t constops_links_score(const uint8_t* req, size_t req_len, ConstopsBuffer* out);

/* Solve a schedule from a `SolveScheduleRequest` JSON payload. */
int32_t constops_schedule_solve(const uint8_t* req, size_t req_len, ConstopsBuffer* out);

/* Validate a schedule from a `ValidateScheduleRequest` JSON payload. */
int32_t constops_schedule_validate(const uint8_t* req, size_t req_len, ConstopsBuffer* out);

/*
 * Ingest an orbit state from an `IngestOrbitStateRequest` JSON payload.
 *
 * Each call uses a fresh ephemeral OrbitRegistry; cross-call lineage is not
 * preserved by the FFI (use the REST API for persistent lineage).
 */
int32_t constops_orbit_ingest(const uint8_t* req, size_t req_len, ConstopsBuffer* out);

/*
 * Release a ConstopsBuffer previously returned by a constops_* call.
 *
 * No-op when `buf.data == NULL && buf.len == 0`. Calling twice on the same
 * non-empty buffer is undefined behaviour.
 */
void constops_buffer_free(ConstopsBuffer buf);

/*
 * Last error message recorded on the calling thread (UTF-8, NUL-terminated).
 *
 * The returned pointer is owned by the FFI and remains valid until the next
 * constops_* call on the same thread or process exit. Returns NULL when no
 * error has been recorded.
 */
const char* constops_last_error(void);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* CONSTOPS_FFI_H */
