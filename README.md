# Reasoning is topologically simple: code and results

This is the repository associated to the paper *Reasoning is topologically simple* (C. Kapulkin and
N. Kershaw). It contains the optimized co-occurrence persistent homology code along with the computational results
of the paper's four experiments.


## Layout

```
cooccurrence_code/
  cooccurrence_ph.jl             co-occurrence persistent homology (Julia)
experimental_results/
  1_reasoning_vs_nonreasoning_prose/   texts.csv, genre_summary.csv, hypothesis_tests.csv
  2_openai_vs_human_papers/            pairs.csv, summary.csv
  3_llm_capability_vs_topology/        models.csv, traces.csv, summary.csv, model_responses/
  4_shade_arena/                       runs.csv, tasks.csv, summary.csv
```

In all results files:
- `total_persistence` is the observed value;
- `null_total_persistence` is the mean over the 200 shuffles;
- `total_persistence_over_null` is their ratio, the quantity used in the paper.

## `cooccurrence_code/`

`cooccurrence_ph.jl` the main function is `cooccurrence_persistence(sets; maxdim=1)`.
- **Input:** a filtration given as a list of sets. Step *t* adds the full simplex on the *t*-th set,
  i.e. the words of the *t*-th sentence. Vertex labels can be words, integers or any other sortable
  values.
- **Output:** the persistence diagram in dimensions 0 to `maxdim`, as sorted `(dimension, birth, death)`
  tuples. Steps are numbered from 1, a death of `-1` marks a class that never dies, and zero-length
  intervals are omitted.

```julia
include("cooccurrence_ph.jl")
using .CooccurrencePH

sentences = [["square", "odd"], ["square", "remainder"],
             ["odd", "remainder"], ["odd", "square", "remainder"]]
cooccurrence_persistence(sentences; maxdim=1)
# [(0, 1, -1), (1, 3, 4)]
```

This is the example in Figure 1 in the paper. 

## `experimental_results/`

### 1. Reasoning vs non-reasoning prose

96 Project Gutenberg texts, 24 per genre (mathematics, argumentative philosophy, history, short
stories). Each text is cut to 300 sentences taken from its middle.

| File | Contents |
|---|---|
| `texts.csv` | one row per text: genre, Gutenberg eBook number, author, title and the three total-persistence values |
| `genre_summary.csv` | median and quartiles per genre (Table 1) |

The texts themselves are not included; each can be found on Project Gutenberg by its eBook number.

### 2. OpenAI vs human mathematics papers

135 pairs, each pairing an OpenAI mathematics paper with a human-written paper on a similar sibject. 

Each paper is reduced to its body: preamble, abstract, bibliography, all mathematics and figures are
removed. Each pair is truncated to the shorter paper's sentence count.

| File | Contents |
|---|---|
| `pairs.csv` | the pairing (OpenAI title; human authors, year, title, venue and arXiv version), the number of sentences compared, both papers' total-persistence values and their difference |
| `summary.csv` | the number and share of pairs with the human paper lower, the one-sided Wilcoxon signed-rank p, and the medians, quartiles, means and standard deviations |

### 3. LLM capability vs chain-of-thought topology

48 open-weight models with an Epoch Capabilities Index (ECI), each run on 30 AIME 2026
problems. Each reasoning trace is truncated to its first 100 sentences, and traces shorter than that
are not scored.

| File | Contents |
|---|---|
| `models.csv` | per model: ECI, median total persistence over null across its traces, AIME accuracy and number of traces scored |
| `traces.csv` | per model and problem: correctness and the three total-persistence values |
| `summary.csv` | Spearman rho between ECI and the model medians (one-sided p), and the within-lab comparison: the share of same-lab model pairs in which the higher-ECI model has the lower median, with its permutation p |
| `model_responses/` | the raw API response for every model and problem: problem, extracted and expected answer, correctness, token usage and the reasoning text |

### 4. SHADE-Arena

26 SHADE-Arena tasks, each run three times in the benign version and three times with the sabotage
side task. Agent `deepseek-v4-flash-0731`. Within each task, every run is truncated to the sentence
count of the task's shortest run. Task names are replaced by hashes.

| File | Contents |
|---|---|
| `runs.csv` | one row per run: task, version (`benign` or `sabotage`), run number, sentence count and the three total-persistence values |
| `tasks.csv` | per task: mean total persistence over null of the benign and of the sabotage runs, and their difference |
| `summary.csv` | the number of tasks in which the benign version is simpler, the overall benign and sabotage means, and the one-sided Wilcoxon signed-rank p |

## Data sources and licences

- `3_llm_capability_vs_topology/model_responses/` contains the AIME 2026 problem statements from the
  MathArena dataset (<https://huggingface.co/datasets/MathArena/aime_2026>). That dataset is licensed
  CC BY-NC-SA 4.0; please cite MathArena when using it.
- The OpenAI papers are from <https://github.com/openai/math>. Human-paper metadata (authors, titles) is from arXiv. No paper text is included.
- The prose texts are from Project Gutenberg (<https://www.gutenberg.org>). Only catalogue metadata is
  included.
- The SHADE-Arena tasks are from J. Kutasov et al., *SHADE-Arena: Evaluating Sabotage and Monitoring in
  LLM Agents* (<https://arxiv.org/abs/2506.15740>).
