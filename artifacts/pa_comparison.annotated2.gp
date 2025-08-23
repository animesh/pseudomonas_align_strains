set terminal pngcairo size 1600,1200 enhanced font "Courier,10"
set output '/home/ash022/pseudomonas_align_strains/pa_comparison.annotated2.png'

set size 1,1
set grid
unset key
set border 10
set tics scale 0
set xlabel "CP011857.1"
set ylabel "QRY"
set format "%.0f"
set xrange [1:6833187]
set yrange [1:7550891]
set style line 1  lt 1 lw 2 pt 6 ps 1
set style line 2  lt 3 lw 2 pt 6 ps 1
set style line 3  lt 2 lw 2 pt 6 ps 1

# Highlight blaVIM2 -> ST matches as horizontal bands (y ranges from bla_ST.coords)
# Match 1: ST 1572..2214
set object 1 rect from 1,1572 to 6833187,2214 fc rgb "red"  fillstyle solid 0.15 noborder
# Match 2: ST 148..1002
set object 2 rect from 1,148  to 6833187,1002 fc rgb "red"  fillstyle solid 0.15 noborder
# Match 3: ST 4197..4478
set object 3 rect from 1,4197 to 6833187,4478 fc rgb "red"  fillstyle solid 0.15 noborder
# Match 4: ST 56..155  (normalized from 155..56)
set object 4 rect from 1,56   to 6833187,155  fc rgb "red"  fillstyle solid 0.15 noborder

# Plot the existing forward/reverse mummerplot data
plot \
 'pa_comparison.fplot' title 'FWD' with lp ls 1, \
 'pa_comparison.rplot' title 'REV' with lp ls 2

# end
