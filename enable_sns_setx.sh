#!/bin/bash

# add 'set -x' after every '#!/bin/bash'

enable_setx () {
    local filepath="$1"
    sed -i -e 's|^\(#!/bin/bash\)$|\1\nset -x|g' "$filepath"
}

scripts_dir="$1"

# check before
find "$scripts_dir" -type f -name "*.sh" -exec grep -A 5 'bin/bash' {} \;

find "$scripts_dir" -type f -name "*.sh" | while read item; do
    echo "$item"
    enable_setx "$item"
done
