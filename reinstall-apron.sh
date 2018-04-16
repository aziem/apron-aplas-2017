#!/bin/bash

opam remove apron -y
opam pin add apron . -y
notify-send -u critical "Done reinstalling APRON"
