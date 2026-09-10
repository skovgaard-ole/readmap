# `read_print_23.pl` — the Perl original

This is the program `readmap` was ported from. It is **not maintained and not
run**; it is kept here as the reference the port was audited against.

`CHANGES.md` documents every way the Python output differs from this program's,
and cites line numbers in this file — for example `read_print_23.pl:1141`, the
line that cut every gene name to its last five characters. Those citations are
only checkable while the file is here.

To read the algorithm rather than the code, `code_structure.md` §5 states it
step by step. `port_plan.md`, beside this file, is the original audit of the
Perl and the plan the Python package was built from — also historical, and also
not maintained.

Requires Perl and Ghostscript, neither of which `readmap` needs.
