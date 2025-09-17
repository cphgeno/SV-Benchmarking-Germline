#!/usr/bin/awk -f
BEGIN { FS=OFS="\t" }

# Only process @SQ lines
$1 == "@SQ" {
    for (i = 1; i <= NF; i++) {
        if ($i ~ /^SN:/) {
            sub(/^SN:/, "", $i)
            gsub(/[*:,\/()]/, "_", $i)
            $i = "SN:" $i
        }
    }
}

# Print all lines (including non-@SQ unchanged)
{ print }

