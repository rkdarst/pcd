
set log x
set style data lines

plot "tmp-polopatribes.txt" \
       using 1:2 t "q",  \
    "" using 1:4 t "H",  \
    "" using 1:5 t "I",  \
    "" using 1:6 t "VI", \
    "" using 1:7 t "In", \
    "nussinov_2009_fig8_H.txt" using 1:2
#    "" using 1:3 t "E",  \

plot [] [] "tmp-polopatribes.txt" \
       using 1:2 lt 1 lw 3 t "q",  \
    "" using 1:4 lt 2 lw 3 t "H",  \
    "" using 1:5 lt 3 lw 3 t "I",  \
    "data/highland_peter/highlandRetestrw.txt" using 4:15 lt 1 lw 1 t "q-peter", \
    "data/highland_peter/highlandRetestrw.txt" using 4:11 lt 2 lw 1 t "H-peter", \
    "data/highland_peter/highlandRetestrw.txt" using 4:9  lt 3 lw 1 t "I-peter"

    "tmp-polopatribes2.txt" \
       using 1:2 lt 4 lw 3 t "q",  \
    "" using 1:4 lt 5 lw 3 t "H",  \
    "" using 1:5 lt 6 lw 3 t "I",  \
    "data/nussinov_2009_fig8_q.txt" using 1:2 lt 1 lw 1 t "q-paper"
    "data/nussinov_2009_fig8_H.txt" using 1:2 lt 2 lw 1 t "H-paper", \
    "data/nussinov_2009_fig8_I.txt" using 1:2 lt 3 lw 1 t "I-paper"

plot [] [-5000:] "tmp-polopatribes.txt" \
       using 1:($3*16) lt 4 lw 3 t "E",  \
    "data/highland_peter/highlandRetestrw.txt" using 4:17 lt 4 lw 1 t "E-peter"


plot "tmp-dolphins.txt" \
       using 1:($2/10) lt 1 lw 3 t "q/10",  \
    "" using 1:4       lt 2 lw 3 t "H",  \
    "" using 1:5       lt 3 lw 3 t "I",  \
    "data/nussinov_2009_fig6_q.txt" using 1:($2/10) lt 1 lw 1 t "(q/10)-paper", \
    "data/nussinov_2009_fig6_H.txt" using 1:2       lt 2 lw 1 t "H-paper", \
    "data/nussinov_2009_fig6_I.txt" using 1:2       lt 3 lw 1 t "I-paper"
#    "" using 1:(($3+160)/20) lt 4 lw 3 t "E",  \


plot "tmp-lattice.txt" \
       using 1:($2/10) lt 1 lw 3 t "q",  \
    "" using 1:(($3+160)/20) lt 4 lw 3 t "E",  \
    "" using 1:4       lt 2 lw 3 t "H",  \
    "" using 1:5       lt 3 lw 3 t "I"

plot "tmp-256node.txt" \
       using 1:($2/10) lt 1 lw 3 t "q/10",  \
    "" using 1:4       lt 2 lw 3 t "H",  \
    "" using 1:5       lt 3 lw 3 t "I",  \
    "data/nussinov_2009_fig2_q.txt" using 1:($2/10) lt 1 lw 1 t "(q/10)-paper", \
    "data/nussinov_2009_fig2_H.txt" using 1:2       lt 2 lw 1 t "H-paper", \
    "data/nussinov_2009_fig2_I.txt" using 1:2       lt 3 lw 1 t "I-paper"
#    "" using 1:(($3+160)/20) lt 4 lw 3 t "E",  \


plot [] [] "tmp-karate.txt" \
       using 1:($2/4) lt 1 lw 3 t "q/4",  \
    "" using 1:4 lt 2 lw 3 t "H",  \
    "" using 1:5 lt 3 lw 3 t "I",  \
    "" using 1:6 lt 4 lw 3 t "VI"
    "data/highland_peter/highlandRetestrw.txt" using 4:15 lt 1 lw 1 t "q-peter", \
    "data/highland_peter/highlandRetestrw.txt" using 4:11 lt 2 lw 1 t "H-peter", \
    "data/highland_peter/highlandRetestrw.txt" using 4:9  lt 3 lw 1 t "I-peter"
set term png size 1024,700 ; set output "karate.png" ; replot ; set term pop