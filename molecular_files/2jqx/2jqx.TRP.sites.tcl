mol new dummy.xyz
topo clearbonds
proc vmd_draw_arrow {mol start end color} {
    # an arrow is made of a cylinder and a cone
    set middle [vecadd $start [vecscale 0.800000 [vecsub $end $start]]]
    graphics $mol color $color
    graphics $mol cylinder $start $middle radius 0.100000
    graphics $mol cone $middle $end radius 0.200000
}
draw arrow {19.413000 26.246000 1.503000} {20.335503 25.940470 1.738882} red
draw arrow {19.413000 26.246000 1.503000} {19.180389 26.293630 2.474403} red
draw arrow {19.413000 26.246000 1.503000} {19.104973 25.295009 1.475869} red
draw arrow {18.156000 28.084000 0.318000} {17.935815 29.047920 0.168414} blue
draw arrow {18.156000 28.084000 0.318000} {17.205834 27.837367 0.127324} blue
draw arrow {18.156000 28.084000 0.318000} {17.935311 28.184147 1.288189} blue
