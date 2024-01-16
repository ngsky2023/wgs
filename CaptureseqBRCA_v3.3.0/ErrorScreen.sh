#!/bin/bash
## run this script under the output directory of NGS processing directory
awk '/error|fail/{print FILENAME, $0}' ./*/*/*.sh.e* > ErrorScreen.txt
