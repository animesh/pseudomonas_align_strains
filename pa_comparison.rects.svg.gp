set terminal svg size 1600,1200 enhanced font "Courier,10"
set output '/home/ash022/pseudomonas_align_strains/pa_comparison.rects.svg'

set size 1,1
set grid
unset key
set border 10
set tics scale 0
set xlabel "CP011857.1"
set ylabel "QRY (ST coordinates)"
set format "%.0f"
set xrange [1:6833187]
set yrange [1:7550891]
set style line 1  lt 1 lw 2 pt 6 ps 1
set style line 2  lt 3 lw 2 pt 6 ps 1
set style line 3  lt 2 lw 2 pt 6 ps 1

# Corrected: highlight blaVIM2 -> ST matches using ST coordinate ranges (S1..E1 from bla_ST.coords)
# Match A: ST 3843..4485  (center ~4164)
set object 1 rect from 1,3843 to 6833187,4485 fc rgb "red"  fillstyle solid 0.15 noborder
# Match B: ST 4477..5332  (center ~4905)
set object 2 rect from 1,4477 to 6833187,5332 fc rgb "red"  fillstyle solid 0.15 noborder
# Match C: ST 1..282      (center ~142)
set object 3 rect from 1,1    to 6833187,282  fc rgb "red"  fillstyle solid 0.15 noborder
# Match D: ST 124..223    (center ~174)
set object 4 rect from 1,124  to 6833187,223  fc rgb "red"  fillstyle solid 0.15 noborder

# Add arrows and labels pointing to each match center (x head at mid-plot)
X_HEAD = (6833187+1)/2.0
X_LABEL = 150000.0

# centers
Y1 = (3843+4485)/2.0
Y2 = (4477+5332)/2.0
Y3 = (1+282)/2.0
Y4 = (124+223)/2.0

set label 1 "blaVIM2" at X_LABEL, Y1+300 tc rgb "blue" font ",12"
set arrow 1 from X_LABEL, Y1+300 to X_HEAD, Y1 nohead lc rgb "blue" lw 2

set label 2 "blaVIM2" at X_LABEL, Y2+300 tc rgb "blue" font ",12"
set arrow 2 from X_LABEL, Y2+300 to X_HEAD, Y2 nohead lc rgb "blue" lw 2

set label 3 "blaVIM2" at X_LABEL, Y3+200 tc rgb "blue" font ",12"
set arrow 3 from X_LABEL, Y3+200 to X_HEAD, Y3 nohead lc rgb "blue" lw 2

set label 4 "blaVIM2" at X_LABEL, Y4-100 tc rgb "blue" font ",12"
set arrow 4 from X_LABEL, Y4-100 to X_HEAD, Y4 nohead lc rgb "blue" lw 2

# Plot the existing forward/reverse mummerplot data
plot \
 'pa_comparison.fplot' title 'FWD' with lp ls 1, \
 'pa_comparison.rplot' title 'REV' with lp ls 2

# end

# Precise ATCC x ST rectangles computed from pa_comparison.coords and bla_ST.coords
# rect 1: ST JAFFXY010000040.1 458264-458906 -> ATCC 1-6833187
set object 101 rect from 1,458264 to 6833187,458906 fc rgb "green" fillstyle solid 0.25 noborder
# rect 2: ST JAFFXY010000040.1 458898-459753 -> ATCC 1-6833187
set object 102 rect from 1,458898 to 6833187,459753 fc rgb "green" fillstyle solid 0.25 noborder
# rect 3: ST JAFFXY010001665.1 342093-342374 -> ATCC 1-6833187
set object 103 rect from 1,342093 to 6833187,342374 fc rgb "green" fillstyle solid 0.25 noborder
# rect 4: ST JAFFXY010001793.1 566688-566787 -> ATCC 1-6833187
set object 104 rect from 1,566688 to 6833187,566787 fc rgb "green" fillstyle solid 0.25 noborder

# Plot the existing forward/reverse mummerplot data
plot \
+ '/home/ash022/pseudomonas_align_strains/pa_comparison.fplot' title 'FWD' with lp ls 1, \
+ '/home/ash022/pseudomonas_align_strains/pa_comparison.rplot' title 'REV' with lp ls 2
