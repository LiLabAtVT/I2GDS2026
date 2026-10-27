```bash
export PATH="$(ls -d ~/.vscode-server/extensions/anthropic.claude-code-*/resources/native-binary | tail -1):$PATH"
claude plugin marketplace add DietrichGebert/ponytail
claude plugin install ponytail@ponytail
claude plugin list
cd /projects/intro2gds/I2GDS2026/Team1/Kamal/Claude_Code
claude -p "/ponytail-help"
claude -p "/ponytail-gain"
claude -p "/ponytail-debt"
```
