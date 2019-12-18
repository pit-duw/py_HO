
  reset
  set lmargin 10
  set rmargin 0
  set terminal postscript col enhanced
  set output "fit_results.ps"
  set style data histeps
  set key off
  set style line 1 lt 1 lc rgb "#006D4F" lw 1.8
  set style line 2 lt 1 lc rgb "#B90091" lw 1.8

  set label "{/Symbol @\326\140}S=   5.02TeV, -0.90<y<  0.90" font ",14" at graph 0.1, graph 0.94
  set title "SET NUMBER =  1 distribution" font "Helvetica, 20"
  set xlabel "P_T [GeV]" font "Helvetica, 20"
  set ylabel "d^2{/Symbol s}/dP_Tdy [nb/GeV]" font "Helvetica, 20"
  set label front "HELAC-ONIA" font "Helvetica, 15" rotate by 90 at graph 1.02, graph 0.4
  set xrange [    0.00000:  10.00000]
  set format y "10^{%T}"
  set key horizontal noreverse maxcols 1 width -4 at graph 0.92, graph 0.9
  set logscale y
  plot \
  "./comparison_1.dat" using 1:($4):($4+$5):($4-$5) with yerror ls 1 title "EX",\
  "./comparison_1.dat" using 1:($2) with histeps ls 2 title "TH"
  unset label

  set label "{/Symbol @\326\140}S=   5.02TeV, -0.91<y<  0.90" font ",14" at graph 0.1, graph 0.94
  set title "SET NUMBER =  2 distribution" font "Helvetica, 20"
  set xlabel "P_T [GeV]" font "Helvetica, 20"
  set ylabel "d^2{/Symbol s}/dP_Tdy [nb/GeV]" font "Helvetica, 20"
  set label front "HELAC-ONIA" font "Helvetica, 15" rotate by 90 at graph 1.02, graph 0.4
  set xrange [    1.30000:  10.00000]
  set format y "10^{%T}"
  set key horizontal noreverse maxcols 1 width -4 at graph 0.92, graph 0.9
  set logscale y
  plot \
  "./comparison_2.dat" using 1:($4):($4+$5):($4-$5) with yerror ls 1 title "EX",\
  "./comparison_2.dat" using 1:($2) with histeps ls 2 title "TH"
  unset label

  set label "{/Symbol @\326\140}S=   8.02TeV, -3.98<y< -2.50" font ",14" at graph 0.1, graph 0.94
  set title "SET NUMBER =  3 distribution" font "Helvetica, 20"
  set xlabel "P_T [GeV]" font "Helvetica, 20"
  set ylabel "d^2{/Symbol s}/dP_Tdy [nb/GeV]" font "Helvetica, 20"
  set label front "HELAC-ONIA" font "Helvetica, 15" rotate by 90 at graph 1.02, graph 0.4
  set xrange [    0.00000:  20.00000]
  set format y "10^{%T}"
  set key horizontal noreverse maxcols 1 width -4 at graph 0.92, graph 0.9
  set logscale y
  plot \
  "./comparison_3.dat" using 1:($4):($4+$5):($4-$5) with yerror ls 1 title "EX",\
  "./comparison_3.dat" using 1:($2) with histeps ls 2 title "TH"
  unset label

  set label "{/Symbol @\326\140}S=   8.02TeV,  1.56<y<  3.06" font ",14" at graph 0.1, graph 0.94
  set title "SET NUMBER =  4 distribution" font "Helvetica, 20"
  set xlabel "P_T [GeV]" font "Helvetica, 20"
  set ylabel "d^2{/Symbol s}/dP_Tdy [nb/GeV]" font "Helvetica, 20"
  set label front "HELAC-ONIA" font "Helvetica, 15" rotate by 90 at graph 1.02, graph 0.4
  set xrange [    0.00000:  20.00000]
  set format y "10^{%T}"
  set key horizontal noreverse maxcols 1 width -4 at graph 0.92, graph 0.9
  set logscale y
  plot \
  "./comparison_4.dat" using 1:($4):($4+$5):($4-$5) with yerror ls 1 title "EX",\
  "./comparison_4.dat" using 1:($2) with histeps ls 2 title "TH"
  unset label

  set label "{/Symbol @\326\140}S=   5.02TeV, -0.90<y<  0.90" font ",14" at graph 0.1, graph 0.94
  set title "SET NUMBER =  5 distribution" font "Helvetica, 20"
  set xlabel "P_T [GeV]" font "Helvetica, 20"
  set ylabel "d^2{/Symbol s}/dP_Tdy [nb/GeV]" font "Helvetica, 20"
  set label front "HELAC-ONIA" font "Helvetica, 15" rotate by 90 at graph 1.02, graph 0.4
  set xrange [    0.00000:  20.10050]
  set format y "10^{%T}"
  set key horizontal noreverse maxcols 1 width -4 at graph 0.92, graph 0.9
  set logscale y
  plot \
  "./comparison_5.dat" using 1:($4):($4+$5):($4-$5) with yerror ls 1 title "EX",\
  "./comparison_5.dat" using 1:($2) with histeps ls 2 title "TH"
  unset label

  set label "{/Symbol @\326\140}S=   5.02TeV, -0.91<y<  0.90" font ",14" at graph 0.1, graph 0.94
  set title "SET NUMBER =  6 distribution" font "Helvetica, 20"
  set xlabel "P_T [GeV]" font "Helvetica, 20"
  set ylabel "d^2{/Symbol s}/dP_Tdy [nb/GeV]" font "Helvetica, 20"
  set label front "HELAC-ONIA" font "Helvetica, 15" rotate by 90 at graph 1.02, graph 0.4
  set xrange [    0.00000:  20.10050]
  set format y "10^{%T}"
  set key horizontal noreverse maxcols 1 width -4 at graph 0.92, graph 0.9
  set logscale y
  plot \
  "./comparison_6.dat" using 1:($4):($4+$5):($4-$5) with yerror ls 1 title "EX",\
  "./comparison_6.dat" using 1:($2) with histeps ls 2 title "TH"
  unset label

  set label "{/Symbol @\326\140}S=   8.02TeV, -3.98<y< -2.50" font ",14" at graph 0.1, graph 0.94
  set title "SET NUMBER =  7 distribution" font "Helvetica, 20"
  set xlabel "P_T [GeV]" font "Helvetica, 20"
  set ylabel "d^2{/Symbol s}/dP_Tdy [nb/GeV]" font "Helvetica, 20"
  set label front "HELAC-ONIA" font "Helvetica, 15" rotate by 90 at graph 1.02, graph 0.4
  set xrange [    0.00000:  20.10050]
  set format y "10^{%T}"
  set key horizontal noreverse maxcols 1 width -4 at graph 0.92, graph 0.9
  set logscale y
  plot \
  "./comparison_7.dat" using 1:($4):($4+$5):($4-$5) with yerror ls 1 title "EX",\
  "./comparison_7.dat" using 1:($2) with histeps ls 2 title "TH"
  unset label

  set label "{/Symbol @\326\140}S=   8.02TeV,  1.56<y<  3.06" font ",14" at graph 0.1, graph 0.94
  set title "SET NUMBER =  8 distribution" font "Helvetica, 20"
  set xlabel "P_T [GeV]" font "Helvetica, 20"
  set ylabel "d^2{/Symbol s}/dP_Tdy [nb/GeV]" font "Helvetica, 20"
  set label front "HELAC-ONIA" font "Helvetica, 15" rotate by 90 at graph 1.02, graph 0.4
  set xrange [    0.00000:  20.10050]
  set format y "10^{%T}"
  set key horizontal noreverse maxcols 1 width -4 at graph 0.92, graph 0.9
  set logscale y
  plot \
  "./comparison_8.dat" using 1:($4):($4+$5):($4-$5) with yerror ls 1 title "EX",\
  "./comparison_8.dat" using 1:($2) with histeps ls 2 title "TH"
