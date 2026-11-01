#!/bin/sh
# dipole_gain_check.sh — 遠方界の利得の検証 (CI 用)
#
# data/sample/dipole_gain.ofd (半波長ダイポール、完全導体) をソルバーと ofd_post で
# 走らせ、far2d.log / far1d.log から次を判定する (根拠と実測値は .ofd のコメント):
#   1. (1/4π)∮G dΩ = 1             (損失が無いので放射電力 = 受け入れた電力)  ±3%
#   2. matchingloss = 1 では 1-|Γ|^2 (Γ は ofd.log の 3 GHz の反射)          ±3%
#   3. 正面 (θ = 90) の利得 2.15 dBi ± 0.4 dB、軸方向 (θ = 0) は -20 dB 未満
#
# 使い方 : dipole_gain_check.sh <ソルバー(絶対パス)> [作業ディレクトリ]
#
# 環境変数 (tpa_slab_check.sh と同じ):
#   OFD_LAUNCHER … ソルバーの前に置くコマンド (例: "mpirun -n 2")
#   OFD_ARGS     … ソルバーに渡す追加オプション (例: "-cpu")
#   OFD_POST     … 後処理 (既定: ソルバーと同じ場所の ofd_post)

set -e

OFD="$1"
WORK="${2:-.}"
LAUNCHER="${OFD_LAUNCHER:-}"
EXTRA="${OFD_ARGS:-}"
SRC="$(cd "$(dirname "$0")" && pwd)/dipole_gain.ofd"

if [ -z "$OFD" ]; then
	echo "Usage: dipole_gain_check.sh <ofd> [workdir]" >&2
	exit 2
fi

POST="${OFD_POST:-$(dirname "$OFD")/ofd_post}"
if [ ! -x "$POST" ] && [ -x "$POST.exe" ]; then
	POST="$POST.exe"
fi
if [ ! -x "$POST" ]; then
	echo "*** ofd_post not found : $POST" >&2
	exit 1
fi

mkdir -p "$WORK/g0" "$WORK/g1"
cp "$SRC" "$WORK/g0/run.ofd"
# $LAUNCHER / $EXTRA は複数語に分解させたいので意図的にクォートしない
(cd "$WORK/g0" && $LAUNCHER "$OFD" $EXTRA -n 2 run.ofd > /dev/null)
grep -q "normal end" "$WORK/g0/ofd.log" || { echo "*** solver failed" >&2; exit 1; }
(cd "$WORK/g0" && "$POST" -n 2 run.ofd > post.txt 2>&1)

# 不整合損を含める : 同じ ofd.out を後処理だけやり直す
awk '/^end/ { print "matchingloss = 1" } { print }' "$SRC" > "$WORK/g1/run.ofd"
cp "$WORK/g0/ofd.out" "$WORK/g1/"
(cd "$WORK/g1" && "$POST" -n 2 run.ofd > post.txt 2>&1)

for d in g0 g1; do
	[ -s "$WORK/$d/far2d.log" ] || { echo "*** $d : far2d.log not written" >&2; exit 1; }
done

# far2d.log : 2 行のヘッダの後に  No. No. theta phi E-abs[dB] ...
integ() {
	awk 'NR > 2 && NF >= 5 {
		th = $3; ph = $4
		if (ph >= 360) next                     # φ = 360 は φ = 0 と同じ点
		w = ((th == 0) || (th == 180)) ? 0.5 : 1 # θ の台形則
		s += w * exp($5 / 10 * log(10)) * sin(th * 3.14159265358979 / 180)
	}
	END { d = 10 * 3.14159265358979 / 180; printf "%.6f\n", s * d * d / (4 * 3.14159265358979) }' "$1"
}
I0=$(integ "$WORK/g0/far2d.log")
I1=$(integ "$WORK/g1/far2d.log")

# 3 GHz の反射 [dB] : インピーダンス表 (Rin の見出しの次から) の 3e9 の行の 6 列目
REF=$(awk '/Rin\[ohm\]/ { t = 1; next } t && ($1 + 0 == 3e9) { print $6; exit }' "$WORK/g0/ofd.log")
[ -n "$REF" ] || { echo "*** reflection at 3 GHz not found in ofd.log" >&2; exit 1; }

# far1d.log (X 面) : 番号 角度 E-abs[dB] ...
G90=$(awk 'NR > 2 && NF > 3 && ($2 + 0 == 90) { print $3; exit }' "$WORK/g0/far1d.log")
G0=$(awk 'NR > 2 && NF > 3 && ($2 + 0 == 0) { print $3; exit }' "$WORK/g0/far1d.log")

awk -v i0="$I0" -v i1="$I1" -v ref="$REF" -v g90="$G90" -v g0="$G0" 'BEGIN {
	status = 0
	t1 = 1 - exp(ref / 10 * log(10))
	ok = ((i0 > 0.97) && (i0 < 1.03))
	printf "radiated / accepted       = %.4f (expected 1)              -> %s\n", i0, ok ? "OK" : "NG"; if (!ok) status = 1
	r = (t1 > 0) ? i1 / t1 : 0
	ok = ((r > 0.97) && (r < 1.03) && (t1 < 0.9))
	printf "with matching loss        = %.4f (expected 1-|G|^2 = %.4f) -> %s\n", i1, t1, ok ? "OK" : "NG"; if (!ok) status = 1
	ok = ((g90 > 2.15 - 0.4) && (g90 < 2.15 + 0.4))
	printf "broadside gain            = %.3f dBi (expected 2.15)        -> %s\n", g90, ok ? "OK" : "NG"; if (!ok) status = 1
	ok = (g0 < -20)
	printf "axial (theta = 0)         = %.1f dBi (expected null)       -> %s\n", g0, ok ? "OK" : "NG"; if (!ok) status = 1
	if (status) print "*** dipole gain check FAILED"; else print "dipole gain check passed"
	exit status
}'
