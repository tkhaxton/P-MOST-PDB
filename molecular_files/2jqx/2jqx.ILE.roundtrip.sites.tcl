mol new 2jqx.ILE.roundtrip.xyz
topo clearbonds
set sel [atomselect top " name A"]
$sel set radius 1.200000
set sel [atomselect top " name B"]
$sel set radius 1.700000
set sel [atomselect top " name D"]
$sel set radius 1.200000
set sel [atomselect top " name E"]
$sel set radius 1.700000
set sel [atomselect top " name F"]
$sel set radius 1.200000
set sel [atomselect top " name G"]
$sel set radius 1.700000
set sel [atomselect top " name I"]
$sel set radius 1.200000
set sel [atomselect top " name J"]
$sel set radius 1.700000
set sel [atomselect top " name K"]
$sel set radius 1.200000
set sel [atomselect top " name L"]
$sel set radius 1.700000
set sel [atomselect top " name M"]
$sel set radius 1.200000
set sel [atomselect top " name P"]
$sel set radius 1.700000
set sel [atomselect top " name Q"]
$sel set radius 1.200000
set sel [atomselect top " name R"]
$sel set radius 1.700000
set sel [atomselect top " name T"]
$sel set radius 1.200000
set sel [atomselect top " name U"]
$sel set radius 1.700000
set sel [atomselect top " name V"]
$sel set radius 1.200000
set sel [atomselect top " name W"]
$sel set radius 1.700000
set sel [atomselect top " name X"]
$sel set radius 1.700000
color Name A red2
color Name B red2
color Name E red3
color Name F red
color Name I blue2
color Name J blue2
color Name M blue
color Name P blue
color Name X black
proc vmd_draw_arrow {mol start end color} {
    # an arrow is made of a cylinder and a cone
    set middle [vecadd $start [vecscale 0.800000 [vecsub $end $start]]]
    graphics $mol color $color
    graphics $mol cylinder $start $middle radius 0.100000
    graphics $mol cone $middle $end radius 0.200000
}
topo addbond 0 19
topo addbond 1 0
topo addbond 2 1
topo addbond 2 20
topo addbond 3 2
topo addbond 4 1
topo addbond 5 4
topo addbond 6 4
topo addbond 7 5
topo addbond 8 0
topo addbond 9 1
topo addbond 10 4
topo addbond 11 5
topo addbond 12 5
topo addbond 13 6
topo addbond 14 6
topo addbond 15 6
topo addbond 16 7
topo addbond 17 7
topo addbond 18 7
draw arrow {15.864000 29.731000 -10.335000} {16.682973 29.481233 -9.818377} red
draw arrow {15.864000 29.731000 -10.335000} {15.784613 28.790038 -10.664070} red
draw arrow {15.864000 29.731000 -10.335000} {16.432314 29.959486 -11.125451} red
draw arrow {15.337000 31.602000 -8.717000} {14.989641 30.770167 -9.149893} blue
draw arrow {15.337000 31.602000 -8.717000} {14.433112 32.021914 -8.798601} blue
draw arrow {15.337000 31.602000 -8.717000} {15.586656 31.964942 -9.614745} blue
