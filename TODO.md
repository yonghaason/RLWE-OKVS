# TODO

## Port the homomorphic-rotation decode onto the current pipeline

**Idea.** Alternative to the current pre-rotated (sequenced band) layout:
the receiver ships two Galois keys (row `+1`, column swap) and the sender
applies `rotate_rows` / `rotate_columns` during the encrypted decode, so the
receiver no longer sends shifted copies. Per the followup note's analysis
(git history: `af2d47c:docs/kkls-followup-note.tex`), the refined output-swap
schedule needs `L/2 + (w-1)/2` receiver rotations under the two keys and cuts
communication by roughly 38% against the pre-rotated layout's U-shaped
optimum at `eps = 0.2`.

**Where the code lives.**
- LOCAL branch `he-rot-version` (bf54c0f, forked from 730f010):
  `rlwe-okvs/rlwe-okvs/rpmt.{cpp,h}` — Galois key transfer in setup, the
  rotation decode loop, `mHalfSlots`/`mNumHalfBatch` bookkeeping. Delta is
  ~500 lines confined to the HE encode/decode stage.
- LOCAL branch `real-okvs-experiments` (095050c): `homdecode_bfv`, a
  standalone pre-rotated vs homomorphic-rotation compare test — useful for
  verifying the port.
- Neither branch is on GitHub anymore; push them before relying on this
  plan from another machine.

**Port plan.** Add a mode flag to `sspmtParams` (pre-rot | hom-rot) and
branch only the setup (Galois keygen/serialize/receive) and the
`encrypted_decode` path in `sspmt.cpp`. The GMW/secret-sharing layer, the
PSO layer, and the benchmark infrastructure are unaffected; `benchmark.sh`
then measures both modes with the same harness.

**Why parked (2026-08-12).** Rotation key-switching consumes noise budget:
the current parameters decrypt with only ~3 bits of margin, and
`he-rot-version` compensates by bumping the last coeff-modulus prime
(50 -> 60 bits). Reviving this needs a modulus re-selection with a security
re-check, then a joint (m/n, span) re-sweep tuned for the rotation mode.
The sweep methodology is in the followup note.
