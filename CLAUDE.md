# CLAUDE.md

## Project
Pipeline for designing RfxCas13d guides against user-specified targets. The pipeline will be reused for future RfxCas13d guide design.

Two modes:
1. **Target design:** given one or more target sequences, design RfxCas13d guides against them.
2. **Proximity mode:** given previously tested PspCas13b guides (with their target positions and results), identify RfxCas13d guides located close to those tested sites on the same target. The proximity window is [N nt, to be confirmed].

Stack: Python (pandas, numpy, biopython) for analysis; plotly for visualisation. Work runs in a Bash/Linux CLI environment. Conda env: [to be confirmed].

## Working agreement
- Before doing any work, repeat my request back to me in one or two sentences so I can confirm you understood it.
- Before executing anything, state what you plan to do, how, and why. Wait for my approval before running commands, writing files, or changing data.
- Keep responses short and concise unless I ask for more.
- If my request is ambiguous, ask one clarifying question before acting.
- If you are unsure, or the available information is not enough to support an answer, say: "I don't have enough information to support this."

## Critical review
- When I ask for a critique of methods or results, check for: assumptions, confounders, off-target risk, guide efficacy assumptions, sample size, and whether the evidence supports the conclusion.
- Say clearly when evidence is weak or a claim is not supported. Do not soften findings to please me.
- Distinguish established results from speculation.

## Lab book
- The lab book is `docs/RfxCas13d-guide-design-workflow.md`.
- Never write to the lab book without a draft review. Always show me the proposed text in chat first, with no file changes, and wait for my edits and approval.
- Scope: the lab book is curated, not a full record. Include only what I decide is important enough. Write at a higher level than the session log: the overall pipeline, its refinements, and its results. Leave out routine steps, dead ends, and minor commands unless I ask for them.
- Session log vs lab book: the session log records everything done in a session. Do not copy session content into the lab book unless I ask.
- Keep entries factual. State the methods (the exact command run, with script, input file, and output file, so I can rerun it), inputs, and results. Mark anything uncertain as uncertain. Do not add interpretation I have not stated or checked.
- Once I approve, write only the approved text, under a dated heading. Do not rewrite earlier content unless I ask.
- When I say "summarise", you may ask me one short question about whether a specific item from this session belongs in the lab book. Do not draft lab book text in this step, and do not write to the file.

## Working directory boundary
- Work only inside `[working directory, to be confirmed]`.
- Never search, list, read, or write files outside this directory, including with `find`, `grep`, `ls`, or `cat`.
- Do not use absolute paths outside this directory, `cd` outside it, or `~` paths.
- If a task seems to need any other file outside the directory, stop and ask me.

## Input data
- Target sequences: provided by me, with paths listed here as they are approved.
- Previously tested PspCas13b guides: provided by me as a table with guide sequence, target coordinates, and activity result. Use only the files I provide.
- Do not guess or infer sequences from memory or from public databases.

## Coding rules
- Python: follow PEP8. Add comments that explain intent, not just what the code does.
- R and Bash: apply equivalent best practices.
- Bash scripts: start with `set -euo pipefail` and include a usage line near the top.
- Prefer vectorised pandas/numpy operations over loops.
- Flag performance concerns for large datasets (memory, runtime, I/O) before running them.

## Outputs
- Scripts in `scripts/`, results in `results/`, figures in `figures/`.
- Never overwrite raw input data. Write derived files to new names.
- Each run writes its parameters (proximity window, thresholds, input file names) alongside its results, so it can be reproduced.

## End-of-session summary
When I say "summarise":
1. Summarise the full session: what I asked for, what was done (files created or changed, commands run, results), and any open items.
2. Give one commit line in the imperative mood, under 72 characters.
3. Append the same entry to `docs/session-log.md`.
4. Do not run `git commit` or `git push` unless I ask.

## Session log
- Log file: `docs/session-log.md`, inside the project directory only.
- At the start of every session, read `docs/session-log.md` before doing anything else. If the file does not exist, create it with this header and tell me:
  `# Session log: RfxCas13d guide design workflow`
  `Entries are appended by the agent at the end of each session. Newest entries are at the bottom.`
- Append only. Never rewrite, reorder, or delete earlier entries.
- Entry format:

  ## YYYY-MM-DD HH:MM (Australia/Melbourne): <short title>
  **Requested:** <what I asked for>
  **Done:** <files created or changed, commands run, results>
  **Open items:** <anything unfinished or needing my decision>
  **Commit:** <commit line, imperative mood, under 72 characters>

- If the log grows past about 300 lines, read only the last three entries at session start and tell me you did so.