```bash
export PATH="$(ls -d ~/.vscode-server/extensions/anthropic.claude-code-*/resources/native-binary | tail -1):$PATH"
cd /projects/intro2gds/I2GDS2026/Team1/Kamal/Claude_Code
cat .claude/settings.json
claude -p "Create a file named blocked.txt containing hello" --disallowedTools "Bash" "Edit" "Write"
```
