set terminal pngcairo size 1600,1200 enhanced font 'Courier,10'
set output '/home/ash022/pseudomonas_align_strains/pa_comparison.rects.png'
set terminal pngcairo size 1600,1200 enhanced font "Courier,10"
set output '/home/ash022/pseudomonas_align_strains/pa_comparison.rects.png'

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
# Match B: ST 4477..5332  (center ~4905)
# Match C: ST 1..282      (center ~142)
# Match D: ST 124..223    (center ~174)
set object 1 rect from 1,3843 to 6833187,4485 fc rgb "orange"  fillstyle solid 0.5 border lc rgb "black" lw 2 front

set object 2 rect from 1,4477 to 6833187,5332 fc rgb "orange"  fillstyle solid 0.5 border lc rgb "black" lw 2 front

set object 3 rect from 1,1    to 6833187,282  fc rgb "orange"  fillstyle solid 0.5 border lc rgb "black" lw 2 front

set object 4 rect from 1,124  to 6833187,223  fc rgb "orange"  fillstyle solid 0.5 border lc rgb "black" lw 2 front

# Add arrows and labels pointing to each match center (x head at mid-plot)
X_HEAD = (6833187+1)/2.0
X_LABEL = 150000.0
DX = 50000.0

# centers
Y1 = (3843+4485)/2.0
Y2 = (4477+5332)/2.0
Y3 = (1+282)/2.0
Y4 = (124+223)/2.0

set label 1 "blaVIM2" at X_LABEL, Y1+300 tc rgb "blue" font ",12"
set arrow 1 from X_LABEL, Y1+300 to (X_LABEL+DX),(Y1+100) nohead lc rgb "blue" lw 2

set label 2 "blaVIM2" at X_LABEL, Y2+300 tc rgb "blue" font ",12"
set arrow 2 from X_LABEL, Y2+300 to (X_LABEL+DX),(Y2+100) nohead lc rgb "blue" lw 2

set label 3 "blaVIM2" at X_LABEL, Y3+200 tc rgb "blue" font ",12"
set arrow 3 from X_LABEL, Y3+200 to (X_LABEL+DX),(Y3+50) nohead lc rgb "blue" lw 2

set label 4 "blaVIM2" at X_LABEL, Y4-100 tc rgb "blue" font ",12"
set arrow 4 from X_LABEL, Y4-100 to (X_LABEL+DX),(Y4-150) nohead lc rgb "blue" lw 2

# Legend (small color swatches)
# orange = fallback ST y-range (no AT match found)
# cyan   = precise mapping from pa_comparison.coords
set object 200 rect from graph 0.02, graph 0.88 to graph 0.06, graph 0.92 fc rgb "orange" fillstyle solid 0.5 border lc rgb "black" lw 2 front
set label 201 "fallback: ST y-range (no AT match)" at graph 0.07, graph 0.90 tc rgb "black" front
set object 202 rect from graph 0.02, graph 0.82 to graph 0.06, graph 0.86 fc rgb "cyan" fillstyle solid 0.6 border lc rgb "black" lw 2 front
set label 203 "mapped via coords" at graph 0.07, graph 0.84 tc rgb "black" front

# Plot the existing forward/reverse mummerplot data
plot \
 'pa_comparison.fplot' title 'FWD' with lp ls 1, \
 'pa_comparison.rplot' title 'REV' with lp ls 2

# end

# Precise ATCC x ST rectangles computed from pa_comparison.coords and bla_ST.coords
# rect 1: ST JAFFXY010000040.1 458264-458906 -> ATCC 1-1
# legend / note
set label 50 "Note: orange = ST global y-range for blaVIM2 (no AT match found); cyan = precise mapping from coords" at graph 0.02, graph 0.96 tc rgb "black" font ",10" front

# rect 1: widen tiny AT width to make visible
set object 101 rect from 1,458264 to 2000,458906 fc rgb "cyan" fillstyle solid 0.6 border lc rgb "black" lw 2 front
# rect 2: ST JAFFXY010000040.1 458898-459753 -> ATCC 1-1
# rect 2: use cyan and widen if degenerate
set object 102 rect from 1,458898 to 2000,459753 fc rgb "cyan" fillstyle solid 0.6 border lc rgb "black" lw 2 front
# rect 3: ST JAFFXY010001665.1 342093-342374 -> ATCC 1-1
# rect 3: use cyan and widen if degenerate
set object 103 rect from 1,342093 to 2000,342374 fc rgb "cyan" fillstyle solid 0.6 border lc rgb "black" lw 2 front
# rect 4: ST JAFFXY010001793.1 566688-566787 -> ATCC 1-1
# rect 4: widen tiny AT width to make visible
set object 104 rect from 1,566688 to 2000,566787 fc rgb "cyan" fillstyle solid 0.6 border lc rgb "black" lw 2 front

# Plot the existing forward/reverse mummerplot data
plot \
 'pa_comparison.fplot' title 'FWD' with lp ls 1, \
 'pa_comparison.rplot' title 'REV' with lp ls 2

# Grey lines marking regions found in ST but not AT (fallback)
set arrow 301 from 1,Y1 to 6833187,Y1 nohead lc rgb "grey" lw 2 dashtype 2 front
set label 301 "present in ST only" at (X_LABEL+DX+2000),(Y1+200) tc rgb "grey" font ",10"

set arrow 302 from 1,Y2 to 6833187,Y2 nohead lc rgb "grey" lw 2 dashtype 2 front
set label 302 "present in ST only" at (X_LABEL+DX+2000),(Y2+200) tc rgb "grey" font ",10"

set arrow 303 from 1,Y3 to 6833187,Y3 nohead lc rgb "grey" lw 2 dashtype 2 front
set label 303 "present in ST only" at (X_LABEL+DX+2000),(Y3+50) tc rgb "grey" font ",10"

set arrow 304 from 1,Y4 to 6833187,Y4 nohead lc rgb "grey" lw 2 dashtype 2 front
set label 304 "present in ST only" at (X_LABEL+DX+2000),(Y4-150) tc rgb "grey" font ",10"
