# Claude Code Token Usage: Tools Setup

Team 1, Virginia Tech · Kamal Sherawat · Manisha Chaudhary
Intro to Genomic Data Science, Fall 2026

We built this project's 16S rRNA pipeline with Claude Code, an AI coding assistant that runs in the terminal. This README covers **my part: the tools we installed to measure, reduce and control how many tokens Claude Code uses.** It explains where each tool lives, how to set it up on the VT ARC cluster, and what output to expect.

Every command in [`tools/`](tools/) was run exactly as written on ARC (node `falcon2`) on 27 Oct 2026, and all of them finished without errors.

---

## The tools at a glance

| Tool | Goal | What it does | Commands |
|---|---|---|---|
| `/context`, `/usage` | Measure | Shows what fills Claude's context window and how much of the plan limit is used | [`tools/context_usage.md`](tools/context_usage.md) |
| **codeburn** | Measure | Turns Claude Code's local session logs into a cost and usage dashboard | [`tools/codeburn.md`](tools/codeburn.md) |
| **graphify** | Shrink | Builds a knowledge graph of the project so Claude can answer questions from the graph instead of reading whole files | [`tools/graphify.md`](tools/graphify.md) |
| **Ponytail** | Shrink | Plugin that makes Claude write the smallest working change and give shorter replies | [`tools/ponytail.md`](tools/ponytail.md) |
| **Permissions** | Protect | Rules that decide what Claude may run or edit without asking | [`tools/permissions.md`](tools/permissions.md) |
| GitDiagram | Shrink | Website that draws an architecture diagram of any GitHub repo: replace `github.com` with `gitdiagram.com` in the URL | none (website) |

**Where to read more:** [`TOKEN_USAGE_GUIDE.md`](TOKEN_USAGE_GUIDE.md) explains each tool in detail. The [`tools/`](tools/) files contain only the commands, ready to paste.

---

## Why this matters

Every time you send Claude a message, it re-sends **the whole conversation**: instructions, every file it has read, every command output. The text is counted in tokens, so a long session costs more with every turn and eventually runs out of room. These tools let you:

1. **Measure** where the tokens go (`/context`, `/usage`, codeburn).
2. **Shrink** what gets sent (graphify, Ponytail, `/clear`, `/compact`).
3. **Protect** your files and other people's while doing it (permissions).

---

## Where the files are

### In this folder

```
Claude_Code/
├── README.md                ← this file
├── TOKEN_USAGE_GUIDE.md     ← detailed guide to every tool
├── tools/                   ← tested, copy-paste commands, one file per tool
│   ├── codeburn.md
│   ├── context_usage.md
│   ├── graphify.md
│   ├── permissions.md
│   └── ponytail.md
├── CLAUDE.md                ← rules Claude reads at the start of every session (written by graphify)
├── .claude/settings.json    ← project hooks: check the graph before searching/reading (written by graphify)
├── .graphifyignore          ← folders graphify skips (results/fastqc/)
└── graphify-out/            ← the knowledge graph
    ├── graph.json           ← the graph itself (what queries read)
    ├── graph.html           ← interactive picture: open it in a browser
    ├── GRAPH_REPORT.md      ← summary: communities, most connected nodes
    ├── manifest.json, cost.json, cache/   ← bookkeeping for incremental updates
    └── 2026-10-27/          ← automatic backup made before the last rebuild
```

### Installed outside this folder (in your home directory)

| What | Location | Created by |
|---|---|---|
| graphify program (Python 3.12 environment) | `~/.venvs/graphify/` | `python3 -m venv` + `pip install graphifyy` |
| graphify skill for Claude Code | `~/.claude/skills/graphify/` | `graphify install` |
| Ponytail plugin (v5.0.0) | `~/.claude/plugins/cache/ponytail/` | `claude plugin install ponytail@ponytail` |
| Plugin and marketplace settings | `~/.claude/settings.json`, `~/.claude/plugins/installed_plugins.json` | `claude plugin ...` |
| Global instructions for all projects | `~/.claude/CLAUDE.md` | graphify |
| codeburn (downloaded on first `npx` run) | `~/.npm/_npx/` | `npx --yes codeburn` |
| Claude Code session logs (what codeburn reads) | `~/.claude/projects/` | Claude Code, automatically |
| `claude` command-line program | `~/.vscode-server/extensions/anthropic.claude-code-*/resources/native-binary/` | the Claude Code VS Code extension |

---

## Getting started on ARC

None of these tools are on the default PATH on ARC. Each file in [`tools/`](tools/) starts with the line that makes its tool available:

| Tool | Setup line |
|---|---|
| graphify | `module load python312` then `source ~/.venvs/graphify/bin/activate` |
| codeburn | `module load nodejs/22.17.1-GCCcore-14.3.0` |
| `claude` (for `/context`, `/usage`, Ponytail, permissions) | `export PATH="$(ls -d ~/.vscode-server/extensions/anthropic.claude-code-*/resources/native-binary \| tail -1):$PATH"` |

Suggested order for a new user:

1. **graphify** ([`tools/graphify.md`](tools/graphify.md)): install it once, then run `graphify update .` after changing code.
2. **Ponytail** ([`tools/ponytail.md`](tools/ponytail.md)): install it once; it stays on in every session.
3. **Permissions** ([`tools/permissions.md`](tools/permissions.md)): set deny rules before letting Claude work in a shared folder.
4. **Measure** with [`tools/context_usage.md`](tools/context_usage.md) and [`tools/codeburn.md`](tools/codeburn.md) whenever a session feels slow or expensive.

---

## What we measured

A snapshot from running the commands on this project (27 Oct 2026):

| Measurement | Result |
|---|---|
| Knowledge graph size | 156 nodes, 154 links, 19 communities, from 17 files |
| graphify benchmark | about **9× fewer tokens** per question than reading the whole project (~9,000 tokens vs ~1,000) |
| Fresh session size (`/context`) | about 27k of a 1M-token window (3%) before any work |
| Monthly spend (`codeburn overview`) | $13.39 so far in October |
| Config health (`codeburn optimize`) | Grade C, with about 25% potential savings. Top suggestions: lower the bash output limit and remove an unused skill |
| Ponytail benchmark (`/ponytail-gain`, the authors' numbers) | 45% fewer output tokens, 26% lower cost |

On codeburn's dashboard, **Exploration** (Claude reading files to understand the code) was our biggest cost. That is the part graphify is meant to reduce.

---

## A token-saving routine

1. Start Claude Code **from the project folder**, not a shared parent folder.
2. Before a new task: `/clear`. To shrink a long task without losing the thread: `/compact`.
3. Ask codebase questions; graphify's hooks make Claude check the graph first.
4. After editing scripts: `graphify update .`
5. Once a week: `npx --yes codeburn optimize` and act on the top suggestion.

---

## Troubleshooting

| Problem | Fix |
|---|---|
| `graphify: command not found` | `module load python312 && source ~/.venvs/graphify/bin/activate` |
| `codeburn: command not found` | It isn't installed globally. Use `npx --yes codeburn ...` after loading the nodejs module |
| `claude: command not found` | Run the `export PATH=...` line above (the extension's version number changes when it updates) |
| `node: command not found` after `module load nodejs \| tail` | Piping `module load` runs it in a subshell, so nothing gets loaded. Run `module load` on its own line |
| graphify hooks fail, or point at `~/.local/bin/graphify` | Activate the environment, then re-run `graphify claude install` in this folder; it rewrites the hook paths |
| graphify says it "refreshed skill" | A newer graphify updates `~/.claude/skills/graphify` automatically. Set `GRAPHIFY_NO_AUTO_REFRESH=1` to prevent that |
| `graphify path` / `explain` finds no node | Names must match exactly. Copy them from `graphify god-nodes` or `GRAPH_REPORT.md` |
| `/permissions isn't available in this environment` | It only works in an interactive session. Start `claude` and type it there, or pass `--disallowedTools` on the command line |
| `codeburn report` never exits | It's an interactive dashboard that refreshes. Press `q` or Ctrl+C |

---

## Safety notes

- **Never use `--dangerously-skip-permissions` on ARC.** Every team's folders sit side by side, and that flag lets Claude edit or delete anything your account can reach.
- **Raw data stays off GitHub.** `raw_data/`, `refs/` and `results/qiime2/` are about 4 GB, with single files up to 1.9 GB, over GitHub's 100 MB limit. They are listed in `.gitignore`, and the scripts can regenerate them.
- **The `claude -p` commands use your Claude subscription.** The `/ponytail-review` and `/ponytail-audit` commands read a lot of code, so run them when you mean to.
- **Only install Ponytail from `DietrichGebert/ponytail`.** That's the author's official source; copies elsewhere may not be safe.

---

## Links

- graphify: `pip install graphifyy` · [app.graphify.com](https://app.graphify.com)
- codeburn: `npx codeburn` (npm)
- Ponytail: [github.com/DietrichGebert/ponytail](https://github.com/DietrichGebert/ponytail)
- GitDiagram: [gitdiagram.com](https://gitdiagram.com)
- Claude Code docs: [docs.claude.com/en/docs/claude-code](https://docs.claude.com/en/docs/claude-code)
