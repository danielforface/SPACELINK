# Research operations and certification v1
## Every pulse
1. Read AGENTS, PROGRAM, OPERATIONS, STATE, SOURCES, HYPOTHESES and active protocols at exact HEAD; inspect open issues, PRs, CI and prior evidence.
2. Declare pulse ID SL-Pnnn, owner/lane, hypothesis/experiment IDs, concrete acceptance and stopping conditions. Use an issue to claim scope before parallel work; avoid collisions.
3. Search primary literature including contrary results. Record query, database, UTC date, filters, inspected hits, DOI/arXiv version, relevant equations/figures, exact access status, license and reason for exclusion. Use citation chaining. Never cite an uninspected title as proof. Record failed retrievals and null searches.
4. Pre-register experiment parameters, units, priors, controls, test oracle, statistical method, tolerance and compute budget before running. Distinguish exploratory from confirmatory analysis; do not tune on held-out data.
5. Implement minimal auditable code; pin dependencies, seeds, environment and input hashes. Preserve raw observations unmodified; store source pointers/checksums rather than unauthorized copies. Derive outputs by script. Record all attempted runs, failures, exclusions, limitations and negative outcomes.
6. Validate dimensional consistency, limiting cases, independent reference, convergence, sensitivity, uncertainty and null controls. A test passing verifies only what it tests. CI verifies reproducibility, not truth. Require explicit review before CERTIFIED or publication claims.
7. Write a full pulse report and update STATE with exact reproducible commands, evidence links, next task, dependencies and blockers. Commit scoped changes, use draft PRs after bootstrap, avoid force pushes and unapproved merges. Never claim execution without evidence.
## Artifact conventions
research/SOURCES.md: stable S-nnn IDs; source records may expand into sources/S-nnn.md.
research/HYPOTHESES.md: H-nnn, causal mechanism, novelty status, falsifier, dependencies, state.
experiments/SL-Ennn/PROTOCOL.md: preregistration; README or RESULTS.md: run summary; code/config/data manifest/results machine-readable; immutable run records under runs/SL-Ennn/run-id/ (large artifacts external, hash + URL).
research/pulses/SL-Pnnn.md: full log; research/decisions/ADR-nnn.md: scientific/engineering decision, alternatives, evidence and reversal condition. Link each claim to source and experiment IDs. Never silently replace prior negative results or alter original data.
## Evidence ladder
HYPOTHESIS -> DERIVED/SIMULATED -> REPRODUCED (independent reference) -> OBSERVED (physical data) -> independently reviewed CERTIFIED (scope-specific). These are distinct evidence types, not interchangeable steps. UNVERIFIED and REJECTED remain visible. Certification requires source provenance, code revision, commands, environment, inputs and hashes, outputs, limitations, reviewer and explicit pass/fail criteria. Novelty is a separate assessment with closest-work matrix, search coverage and residual uncertainty. No claims of exhaustive coverage or world-first based on search absence.
## Reproducibility and ethics
SI units by default, explicit unit conversions and conventions for particle energies. Public data only with license/terms review. Do not publish secrets, private records, restricted technology, hazardous apparatus instructions or third-party copyrighted text. Human physical experimentation only under appropriate qualified oversight. Report conflicts, AI assistance and uncertainties. Owner decides license and publication/affiliation policy.

## Mandatory pulse closeout
Every pulse report starts with a concise outcome ledger: achieved, disproved/negative, failed/blocked, evidence/CI, scientific claim limits and next task. Deliver the same substance to the owner in the final response. This is a standing owner requirement, not optional.

## Discovery gate (owner direction 2026-11-25)
Before implementation, identify the desired observable, causal chain, new or altered assumption versus prior work, physical constraints, closest known work and cheapest falsifying test. ADR-003 governs ordering. For exploratory pulses, a defensible quantitative elimination or cross-disciplinary hypothesis map is a valid deliverable; code is not mandatory. Record open questions without implying success. Pulse closeout includes a concise big-picture verdict: new feasible mechanism / new testable hypothesis / scoped exclusion / no movement, with evidence and caveats.
