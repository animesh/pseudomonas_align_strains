set output 'pa_comparison.pdf'
set terminal pdfcairo size 11in,8in enhanced font 'Arial,12'
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
set xrange [1:6833187]
set yrange [1:*]
set style line 1  lt 1 lw 2 pt 6 ps 1
set style line 2  lt 3 lw 2 pt 6 ps 1
set style line 3  lt 2 lw 2 pt 6 ps 1
plot \
 "pa_comparison.fplot" title "FWD" w lp ls 1, \
 "pa_comparison.rplot" title "REV" w lp ls 2

