# Extract header
grep "^[@#]" rama.xvg > header.xvg

# Split into chains (here you specify number of rows containing data for each chain, i.e., number of residues)
awk -v Ares=436 -v Bres=427 '
$0 !~ /^[@#]/ {
    count++;
    if (count <= Ares) {
        print $0 >> "rama_chainA_data.xvg";
    } else if (count <= Ares+Bres) {
        print $0 >> "rama_chainB_data.xvg";
    }
    if (count == Ares+Bres) {
        count=0; # reset at end of frame
    }
}
' rama.xvg

# Add header back
cat header.xvg rama_chainA_data.xvg > rama_chainA.xvg
cat header.xvg rama_chainB_data.xvg > rama_chainB.xvg
