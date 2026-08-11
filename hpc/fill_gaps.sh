#!/bin/bash
# Monash Myeloma Model - move spare bootstrap replicates into the failed indices.
#
# Some replicates fail at the risk-equation stage (thin cells in an unlucky resample). Spares are
# generated at 501+ rather than re-seeded in place, because risk_equations.do is deterministic given
# its MI file - the same index would fail the same way. This renames the spares down into the gaps so
# the set is a contiguous 1..500 again.
#
# Replicate b is ONE resample. Its MI file, its coefficients and its three transport artefacts must
# move TOGETHER, or b would pair a BCR prediction from one resample with parameters from another -
# silently, since the filenames would still look right. That is why this moves five files per gap and
# verifies afterwards.
#
# RUN THIS BEFORE THE SCENARIO SIMS. Once sims exist they are a sixth artefact and must move too.
#
# Usage:  bash hpc/fill_gaps.sh            # dry run, shows the mapping, changes nothing
#         bash hpc/fill_gaps.sh --apply    # performs the moves
set -u

ROOT=$HOME/em76/adam
DATA=251128
ANALYSIS=transport_dvd

COEF="$ROOT/analyses/$ANALYSIS/coefficients/bootstrap"
ATR="$ROOT/analyses/$ANALYSIS/outcomes/A_trial/bootstrap"
BTR="$ROOT/analyses/$ANALYSIS/outcomes/B_transport/bootstrap"
MI="$ROOT/data/$DATA/bootstrap"

APPLY=0
[ "${1:-}" = "--apply" ] && APPLY=1

# Gaps in 1..500, and spares available above 500, both keyed on the COEFFICIENTS - that is the
# artefact that actually failed, and the one the sims need.
GAPS=$(for i in $(seq 1 500);   do [ -f "$COEF/coefficients_${ANALYSIS}_B$i.mmat" ] || echo "$i"; done)
SPARES=$(for i in $(seq 501 600); do [ -f "$COEF/coefficients_${ANALYSIS}_B$i.mmat" ] && echo "$i"; done)

NG=$(echo $GAPS   | wc -w)
NS=$(echo $SPARES | wc -w)
echo "gaps in 1..500 : $NG  ->  $GAPS"
echo "spares  > 500  : $NS  ->  $SPARES"
echo

if [ "$NS" -lt "$NG" ]; then
	echo "NOT ENOUGH SPARES: $NS available for $NG gaps. Generate more before applying." >&2
	exit 1
fi

move_one() {   # $1 = from index, $2 = to index
	local s=$1 g=$2 f t
	for pair in \
		"$MI/MRDR Long MI B$s.dta|$MI/MRDR Long MI B$g.dta" \
		"$COEF/coefficients_${ANALYSIS}_B$s.mmat|$COEF/coefficients_${ANALYSIS}_B$g.mmat" \
		"$ATR/bcr_vd_l2_B$s.mmat|$ATR/bcr_vd_l2_B$g.mmat" \
		"$ATR/bcr_dvd_l2_B$s.mmat|$ATR/bcr_dvd_l2_B$g.mmat" \
		"$BTR/transport_dvd_B$s.mmat|$BTR/transport_dvd_B$g.mmat"
	do
		f="${pair%%|*}"; t="${pair##*|}"
		if [ ! -f "$f" ]; then echo "    MISSING SOURCE: $f" >&2; return 1; fi
		if [ "$APPLY" = "1" ]; then mv -f "$f" "$t"; fi
	done
	return 0
}

set -- $SPARES
USED=""
for g in $GAPS; do
	s=${1:-}
	[ -z "$s" ] && { echo "ran out of spares at gap $g" >&2; break; }
	shift
	echo "  B$s -> B$g  (MI, coefficients, bcr_vd, bcr_dvd, transport)"
	move_one "$s" "$g" || { echo "  ABORTED on B$s" >&2; exit 1; }
	USED="$USED $s"
done
echo

# Any spare not used is deleted, so a later glob cannot pick up an orphan replicate that has
# coefficients but no matching transport artefact.
# Exclude the ones just consumed: under --apply they are already gone, but on a DRY RUN they are
# still present and would otherwise be listed as deletions when they are actually moves.
LEFT=""
for i in $(seq 501 600); do
	[ -f "$COEF/coefficients_${ANALYSIS}_B$i.mmat" ] || continue
	case " $USED " in *" $i "*) continue ;; esac
	LEFT="$LEFT $i"
done
if [ -n "$LEFT" ]; then
	echo "unused spares to remove: $LEFT"
	if [ "$APPLY" = "1" ]; then
		for i in $LEFT; do
			rm -f "$MI/MRDR Long MI B$i.dta" \
			      "$COEF/coefficients_${ANALYSIS}_B$i.mmat" \
			      "$ATR/bcr_vd_l2_B$i.mmat" "$ATR/bcr_dvd_l2_B$i.mmat" \
			      "$BTR/transport_dvd_B$i.mmat"
		done
	fi
fi
echo

# VERIFY: every index 1..500 must have all four artefacts the sims need. Silence is the pass.
echo "verifying 1..500 ..."
BAD=0
for i in $(seq 1 500); do
	for f in "$COEF/coefficients_${ANALYSIS}_B$i.mmat" \
	         "$ATR/bcr_vd_l2_B$i.mmat" "$ATR/bcr_dvd_l2_B$i.mmat" \
	         "$BTR/transport_dvd_B$i.mmat"; do
		[ -f "$f" ] || { echo "  B$i missing $(basename "$f")"; BAD=$((BAD+1)); }
	done
done
if [ "$APPLY" = "1" ]; then
	[ "$BAD" -eq 0 ] && echo "  PASS - 500 complete replicates, contiguous." || echo "  $BAD missing artefacts."
else
	echo "  (dry run - rerun with --apply, then this check is the one that matters)"
fi
