A reference simulation can be found in Read MFBO_NMPC\single_simu_NMPC.m

## Working with me

- If a request likely contradicts my underlying intent, or a literal reading would produce something I probably don't want, ask me first instead of implementing it. Surface the conflict and the trade-off; do not silently pick an interpretation.
- To suppress a redundant BibTeX field (e.g. a `url` that duplicates a `doi`), rename it to `xurl` rather than commenting it with `%`. Classic BibTeX does not support `%` comments inside an entry and fails to compile.

## MATLAB style

- Use `%% Section title` for section headers (MATLAB code-folding sections).
- Use `% comment` for inline / explanatory comments.
- Do NOT use decorative banner comments such as `% --- Foo ---------------------` or `% =====================`.



## Academic writing style

Voice: Maintain an objective, authoritative, and precise tone appropriate for peer-reviewed journals. Avoid hyperbolic language, rhetorical questions, anthropomorphism of models/controllers/algorithms (e.g., "the model sees/knows/wants"), and AI-isms (e.g., "delve," "testament," "tapestry," "this is not X, it's Y" / "X is not just Y" negation framing, and labeled mini-headers such as `\paragraph{Short Label.}` or `\emph{Term (parenthetical gloss)}` used to fake structure instead of writing a normal topic sentence).
Structure: Ensure paragraphs have a clear topic sentence, supporting evidence, and a concluding transition.Formatting: Output only in clean Markdown. Format all in-text citations in APA style (or your specified citation style) and ensure they map directly to my sources.Iterative Process: Do not rewrite whole sections at once unless requested. Suggest changes and wait for my approval before finalizing.

- Do not use the "X rather than Y" / "not X, but Y" construction that introduces a frame the reader never held just to negate it (e.g. "we treat these as admissible parametrisations rather than a test of dimensionality reduction"). State the positive claim directly, or use "we do A because of B" / "due to B, we choose A" and lead toward A from the start.
- Do not describe things with loaded terms that assert a conclusion the work is itself examining (e.g. "weakly identified weights"). Describe them literally by what they are (e.g. "the input-magnitude penalty is removed and the volume state weight is fixed").
- Edit surgically. When only one sentence carries the problem, change that sentence and leave the surrounding prose verbatim; do not reword passages that are already precise and in my voice.

### Sentence quality: what makes a bad sentence

A strong sentence carries new information, names its agent, and lets evidence imply the verdict. A weak one carries attitude, a drumroll, or a bare pointer instead. Watch for these faults, ranked by how much they matter to me:

1. **Superlatives / editorial verdicts (worst).** Do not announce a judgment the reader should infer from the evidence. Avoid "the strongest argument," "poorly suited," "unavoidably." State the mechanism and let the verdict follow.
   - Bad: "The strongest argument for the NMPC approach over PID lies in how each responds to..."
   - Bad: "PID control was considered for the single-culture biomass problem but is poorly suited to it."

2. **Vague passives & abstract noun stacks.** Name the agent. Avoid bloodless abstractions that name a thing without saying what it is for.
   - Bad: "Second, the application need was not perceived." (perceived by whom?)

3. **Rhetorical flourish / analogy out of register.** Say it in literal terms, not metaphor or drama. Same point, plain words. This includes single-word choices: prefer a verb whose only meaning is the technical one over a verb that also carries a physical/metaphorical sense.
   - Bad: "Both barriers have now fallen." → literal: "Two recent measurement advances now remove this limitation."
   - Bad: "cutting control variation by a factor of..." (cutting implies slicing with an instrument) → literal: "control variation ... times lower" (no meaning besides the quantitative relation).

4. **Anthropomorphism.** Do not give models, controllers, filters, or algorithms perception, knowledge, or intent. State what the system does mechanically.
   - Bad: "the in- and out-flows the model sees" → "the in- and out-flows applied to the model". Also avoid "the model knows / wants / believes / tries to".

5. **Jargon pile-ups.** Do not stack four or more nouns into a compound; it reads like a variable name. Unstack into a clause.
   - Bad: "the robustness-to-hydrolysate-variability property the project requires."

6. **Bluntness that is actually false.** Do not overstate; a careful reader catches it and it costs credibility.
   - Bad: "the measurements did not exist" (they existed but were hard to obtain and infrequent).

7. **Circularity / redundancy (milder).** Do not restate a premise or a project goal as if the preceding sentences had proven it.

8. **Deictic / throat-clearing openers (mild, avoidable).** "This is the...", "This confirms that...". Not terrible, but prefer to cut the pointer and state the content directly.

**Overclaiming exception — keep the strong conclusion when reporting actual experimental results.** Verbs like *confirm*, *demonstrate*, *prove* must be earned by the evidence cited in the same passage. Describing a mechanism (e.g. an EKF's measurement-switching structure) does not "confirm" achievability, so use "enables" there; experimental results that measured the outcome *do* confirm it, so keep the strong verb.

Also: remove em-dashes (`---`) from documents, replacing with a comma, colon, or semicolon as the context requires.
