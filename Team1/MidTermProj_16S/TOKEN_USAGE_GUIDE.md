# Claude Code Token Usage: Tools and Methods

Team 1, Virginia Tech · Kamal Sherawat · Manisha Chaudhary
Intro to Genomic Data Science, Fall 2026

This guide covers the tools from our "Claude Code Token Usage" presentation that we used while building this 16S rRNA pipeline (`scripts/`, `metadata/`, `results/`). For each one it explains what the tool does, how to install it, and how to use it.

## Why token usage matters

Claude Code works by sending the whole conversation to the model on every turn: the system prompt, tool definitions, memory files (`CLAUDE.md`), every message, and every file or command output Claude has read. That text is split into **tokens**, and all of it is sent again with each request. A long session that has read many files gets slower and more expensive with every turn, and eventually fills the context window.

The tools below fall into three groups, matching the last slide:

| Goal | Tools |
|---|---|
| **Measure:** see where tokens go | `/context`, `/usage`, `codeburn` |
| **Shrink:** send fewer tokens | `/clear`, `/compact`, `graphify`, `ponytail`, GitDiagram |
| **Protect:** limit what Claude can touch | `/permissions`, `.claude/settings.json` |

---

## 1. Measure

### 1.1 `/context`: what fills the context window

**Install:** nothing. It is built into Claude Code.

**Use:** type this in the Claude Code prompt:

```
/context
```

It shows the model, total tokens used out of the window (for example `106.6k / 1.0M tokens (11%)`), and a breakdown by category:

| Category | What it is | How to reduce it |
|---|---|---|
| System prompt / System tools | Claude Code's built-in instructions and tool definitions | Fixed cost, but disabling unused MCP servers removes their tools |
| MCP server instructions | Text added by connected MCP servers | Disconnect servers you don't need (`/mcp`) |
| Memory files | `CLAUDE.md` files (project and `~/.claude/CLAUDE.md`) | Keep them short and specific |
| Skills | Descriptions of installed skills and plugins | Uninstall skills you don't use |
| Messages | The conversation so far, including file contents and command output | `/clear`, `/compact`, read less (see graphify) |
| Autocompact buffer | Space reserved so Claude can summarize before the window fills | Not reducible |

The **Messages** row is the one that grows. Check `/context` whenever a session feels slow.

### 1.2 `/usage`: plan limits and what drives them

**Install:** nothing. It is built into Claude Code.

**Use:**

```
/usage
```

This shows how much of your plan's usage limits you have consumed and what is driving it. Use it to decide when to start a fresh session.

### 1.3 `codeburn`: your session logs as a dashboard

`codeburn` reads Claude Code's local session logs (in `~/.claude/projects/`) and turns them into a terminal dashboard. Nothing is sent anywhere.

**Install:** requires Node.js (which provides `npx`). On a cluster without Node, load it first (for example `module load nodejs`, depending on what your system provides), or run codeburn on your laptop.

```bash
npx codeburn          # run without a global install
```

**Use:**

```bash
codeburn overview     # the dashboard (use `npx codeburn overview` if not installed globally)
```

The dashboard panels:

- **Daily Activity:** cost and number of calls per day.
- **By Project:** cost per project folder (this project shows as `…/Claude_Code`).
- **By Activity:** cost split into Exploration, Coding, and Build/Deploy turns.
- **By Model:** cost, cache-hit rate, and calls per model.
- **Core Tools / Shell Commands:** which tools and shell commands Claude ran most (here `Bash`, `head`, `tail`, `ls`).
- **MCP Servers, Skills & Agents, Workflow:** usage of extensions and corrections/rework rates.

For this project, Exploration (Claude reading files to understand the pipeline) was the biggest cost. That is what graphify reduces.

---

## 2. Shrink

### 2.1 `/clear` and `/compact`: start fresh or summarize

**Install:** nothing. Both are built into Claude Code.

```
/clear       # drop the whole conversation and start with an empty context
/compact     # replace the conversation with a summary, then continue
/compact keep the QIIME2 parameters and the list of failing samples
             # optional: tell it what the summary must keep
```

When to use which:

- **`/clear`** when the task changes (for example, from debugging `03_qiime2.sh` to writing figures in `04_stats.R`). The old conversation is no longer useful and would be resent on every turn.
- **`/compact`** when you are midway through a task and still need the gist of what was done.

The rule from the last slide: **start fresh when the task changes.**

### 2.2 graphify: a knowledge graph of the project

Without a map, Claude answers questions about a project by grepping and reading whole files, and all of that text stays in the context. graphify builds a **knowledge graph** of the project once (files → extract → graph → cluster → report). Claude then queries the graph and gets back only the small relevant subgraph, then reads just the lines it needs. The presentation cites up to a 71% token reduction.

This project's graph is in [`graphify-out/`](graphify-out/):

- `graph.json`: the graph (135 nodes, 139 edges, 13 communities)
- `graph.html`: interactive visualization, open it in a browser
- `GRAPH_REPORT.md`: summary with communities (for example "QIIME2 Processing Script", "Stats & ML Analysis (04_stats.R)", "Co-occurrence Network (05_network.R)") and the most connected "god nodes"

**Install:**

```bash
pip install graphifyy        # note the double "y": the PyPI package name
graphify install             # installs the /graphify skill into Claude Code
graphify claude install      # adds graphify rules to this project's CLAUDE.md and hooks
```

`graphify claude install` is what created this project's [`CLAUDE.md`](CLAUDE.md) rules and the hooks in [`.claude/settings.json`](.claude/settings.json). Those hooks run `graphify hook-guard` before Claude's `Bash`, `Grep`, `Read` and `Glob` calls, reminding it to query the graph before reading raw files.

graphify needs no API key. Code is parsed locally (AST). Documents such as `paper.txt` are summarized by the host agent itself.

**Build the graph** (inside Claude Code, from the project folder):

```
/graphify .
```

**Exclude folders** by listing them in `.graphifyignore` (same syntax as `.gitignore`). This project excludes the 48 FastQC reports:

```
results/fastqc/
```

**Keep it current** after changing code (AST only, no LLM cost):

```bash
graphify update .
```

**Ask questions** (from Claude Code or the terminal):

```bash
graphify query "How are reads denoised in QIIME2?"     # relevant subgraph for a question
graphify query "..." --budget 1500                     # cap the answer at ~1500 tokens
graphify path "03_qiime2.sh" "04_stats.R"              # how two parts connect
graphify explain "SparCC co-occurrence network"        # plain-language explanation of one node
```

Useful `/graphify` options: `--update` (re-extract only changed files), `--wiki` (one markdown page per community), `--svg` / `--graphml` (exports for GitHub or Gephi).

### 2.3 GitDiagram: an architecture diagram from a GitHub repo

GitDiagram ([gitdiagram.com](https://gitdiagram.com)) generates an interactive architecture diagram from any public GitHub repository. It is a website, so there is nothing to install.

**Use:** in the repository URL, replace `github.com` with `gitdiagram.com`:

```
https://github.com/LiLabAtVT/I2GDS2026
https://gitdiagram.com/LiLabAtVT/I2GDS2026
```

The diagram shows components (scripts, modules) and how they connect, and you can zoom and export it. Use it to get the overall picture of an unfamiliar repo before asking Claude about it. You can then point Claude at the specific files you care about instead of having it explore everything. Private repos need a GitHub token entered on the site.

### 2.4 Ponytail: shorter code and shorter replies

Ponytail is a Claude Code plugin ("lazy senior dev mode") that makes Claude write the smallest change that solves the task and reply concisely. Less generated code and text means fewer output tokens, and fewer tokens resent on later turns. The authors' benchmark reports about 45% fewer tokens.

**Install** (in Claude Code, as two separate prompts):

```
/plugin marketplace add DietrichGebert/ponytail
```
```
/plugin install ponytail@ponytail
```

Only install it from `DietrichGebert/ponytail` on GitHub. We use version 5.0.0.

**Use:** once installed it is active every session. Commands:

| Command | What it does |
|---|---|
| `/ponytail lite \| full \| ultra \| off` | Set how aggressive it is (`full` is the default) or turn it off |
| `/ponytail-review` | Review the current change for bugs, security, missing tests, speed, and what to cut |
| `/ponytail-audit` | The same check for the whole repo, most important first |
| `/ponytail-debt` | List the `ponytail:` shortcut comments it deferred |
| `/ponytail-help` | Quick reference |

The last slide also mentions **Headroom** (compresses unnecessary lines) as another plugin in this space. We did not install it for this project.

---

## 3. Protect

### 3.1 Permission rules: `/permissions`

Permission rules control which tools Claude may use, and on what, without asking you first.

**Install:** nothing. It is built into Claude Code.

**Use:**

```
/permissions
```

This opens the rules panel with three lists:

- **Allow:** runs without asking (for example `Bash(ls *)`, `Read`).
- **Ask:** always asks first.
- **Deny:** never allowed. **Deny rules always override allow rules.**

A rule is a tool name, optionally followed by a specifier in parentheses, for example `WebFetch`, `Bash(rm *)`, `Edit(/projects/other_team/**)`. When adding a rule you choose where it is saved: this project just for you (`.claude/settings.local.json`), this project for everyone (`.claude/settings.json`), or all your projects (`~/.claude/settings.json`). Choosing "Always allow" on a permission prompt also saves a rule here.

The same rules can be written directly in a settings file:

```json
{
  "permissions": {
    "allow": ["Bash(ls *)", "Bash(sbatch *)"],
    "deny":  ["Bash(rm -rf *)", "Edit(/projects/intro2gds/I2GDS2026/Team2/**)"]
  }
}
```

### 3.2 Block edits outside your project

Claude Code's **workspace** is the folder you started it in (shown in `/permissions` under Workspace). It can only edit files outside that folder if you approve it. On a shared cluster like ours, where every team's folders sit next to each other, keep it that way:

- Start Claude Code from your own project folder, not from a shared parent.
- Add deny rules for other teams' folders and for destructive commands.
- **Avoid `--dangerously-skip-permissions`.** It turns off every permission prompt, so Claude can edit or delete anything your account can reach, including shared data.

---

## Quick reference

```
MEASURE    /context    /usage    npx codeburn overview
SHRINK     /clear      /compact  /graphify .   graphify update .   graphify query "..."
           gitdiagram.com/<owner>/<repo>       /ponytail [lite|full|ultra]
PROTECT    /permissions          (never --dangerously-skip-permissions)
```
