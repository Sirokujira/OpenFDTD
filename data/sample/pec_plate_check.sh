#!/bin/sh
# pec_plate_check.sh — 平面波の全電界 (散乱界 + 入射界) の検証 (CI 用)
#
# data/sample/pec_plate.ofd (周期境界で無限に広い完全導体の平板に平面波を垂直入射)
# をソルバーと ofd_post で走らせ、near1d.log (導体の法線上の線) と near2d.log
# (法線を含む面) を解析解と比べる
# (根拠と実測値は .ofd のコメント):
#   全電界 (near1d/2dnoinc = 0) : z > 0 で |E| = 2|sin kz| (±0.1)、導体の中で 0 (< 1e-3)
#   散乱界 (near1d/2dnoinc = 1) : 全域で |Es| = 1 (±0.03)
#
# 使い方 : pec_plate_check.sh <ソルバー(絶対パス)> [作業ディレクトリ]
#
# 環境変数 (dipole_gain_check.sh と同じ):
#   OFD_LAUNCHER … ソルバーの前に置くコマンド (例: "mpirun -n 2")
#   OFD_ARGS     … ソルバーに渡す追加オプション (例: "-cpu")
#   OFD_POST     … 後処理 (既定: ソルバーと同じ場所の ofd_post)

set -e

OFD="$1"
WORK="${2:-.}"
LAUNCHER="${OFD_LAUNCHER:-}"
EXTRA="${OFD_ARGS:-}"
SRC="$(cd "$(dirname "$0")" && pwd)/pec_plate.ofd"

if [ -z "$OFD" ]; then
	echo "Usage: pec_plate_check.sh <ofd> [workdir]" >&2
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

mkdir -p "$WORK/total" "$WORK/scat"
cp "$SRC" "$WORK/total/run.ofd"
# $LAUNCHER / $EXTRA は複数語に分解させたいので意図的にクォートしない
(cd "$WORK/total" && $LAUNCHER "$OFD" $EXTRA -n 2 run.ofd > /dev/null)
grep -q "normal end" "$WORK/total/ofd.log" || { echo "*** solver failed" >&2; exit 1; }
(cd "$WORK/total" && "$POST" -n 2 run.ofd > post.txt 2>&1)

# 散乱界 : 同じ ofd.out を near1dnoinc = near2dnoinc = 1 で後処理だけやり直す
sed -e 's/^near1dnoinc = .*/near1dnoinc = 1/' -e 's/^near2dnoinc = .*/near2dnoinc = 1/' "$SRC" > "$WORK/scat/run.ofd"
cp "$WORK/total/ofd.out" "$WORK/scat/"
(cd "$WORK/scat" && "$POST" -n 2 run.ofd > post.txt 2>&1)

for d in total scat; do
	for f in near1d near2d; do
		[ -s "$WORK/$d/$f.log" ] || { echo "*** $d : $f.log not written" >&2; exit 1; }
	done
done

# near1d.log : 2 行のヘッダの後に  No. X Y Z |E| Ex ...      (Z は 4 列目、|E| は 5 列目)
# near2d.log : 2 行のヘッダの後に  No. No. X Y Z |E| Ex ...  (Z は 5 列目、|E| は 6 列目)
# 同じ行数の全電界と散乱界のログを横に並べ、散乱界の |E| は全電界の列数 + |E| の列
check() {
	paste "$WORK/total/$1.log" "$WORK/scat/$1.log" | awk -v name="$1" -v zc="$2" -v ec="$3" '
	NR > 2 && NF > 16 {
		h = NF / 2
		z = $zc; e = $ec; es = $(h + ec)
		if (z > 1e-6) {
			t = 2 * sin(2 * 3.14159265358979 * z / 0.1); if (t < 0) t = -t
			d = e - t; if (d < 0) d = -d
			if (d > front) front = d
			nf++
		}
		else if (z < -1e-6) {
			if (e > inside) inside = e
			ni++
		}
		d = es - 1; if (d < 0) d = -d
		if (d > scat) scat = d
	}
	END {
		status = 0
		if ((nf < 30) || (ni < 15)) { printf "*** %s : too few points (front %d, inside %d)\n", name, nf, ni; status = 1 }
		ok = (front < 0.1)
		printf "%s total, z > 0  : max| |E| - 2|sin kz| | = %.4f (< 0.1)   -> %s\n", name, front, ok ? "OK" : "NG"; if (!ok) status = 1
		ok = (inside < 1e-3)
		printf "%s total, in PEC : max|E| = %.2e (< 1e-3)            -> %s\n", name, inside, ok ? "OK" : "NG"; if (!ok) status = 1
		ok = (scat < 0.03)
		printf "%s scattered     : max| |Es| - 1 | = %.4f (< 0.03)        -> %s\n", name, scat, ok ? "OK" : "NG"; if (!ok) status = 1
		exit status
	}'
}

status=0
check near1d 4 5 || status=1
check near2d 5 6 || status=1
if [ "$status" -ne 0 ]; then
	echo "*** PEC plate check FAILED"
else
	echo "PEC plate check passed"
fi
exit $status
