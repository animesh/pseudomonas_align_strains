# Auto-inserted header: non-interactive terminal + output for headless environments
set terminal pngcairo size 1600,1200 enhanced font "Courier,12"
set output 'ATCC_vs_ST_plot.png'
set ytics ( \
 "" 0 \
)
set size 1,1
set grid
unset key
set border 10
set tics scale 0
set xlabel "CP011857.1"
set ylabel "QRY"
set format "%.0f"
 # removed interactive mouse settings for headless rendering
set xrange [1:6833187]
set yrange [1:*]
set style line 1  lt 1 lw 3 pt 6 ps 1
set style line 2  lt 3 lw 3 pt 6 ps 1
set style line 3  lt 2 lw 3 pt 6 ps 1
plot \
 "ATCC_vs_ST_plot.fplot" title "FWD" w lp ls 1, \
 "ATCC_vs_ST_plot.rplot" title "REV" w lp ls 2
