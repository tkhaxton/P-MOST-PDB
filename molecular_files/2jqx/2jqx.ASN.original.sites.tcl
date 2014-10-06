mol new 2jqx.ASN.original.xyz
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
color Name G red
color Name L blue3
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
topo addbond 0 14
topo addbond 1 0
topo addbond 2 1
topo addbond 2 15
topo addbond 3 2
topo addbond 4 1
topo addbond 5 4
topo addbond 6 5
topo addbond 7 5
topo addbond 8 0
topo addbond 9 1
topo addbond 10 4
topo addbond 11 4
topo addbond 12 7
topo addbond 13 7
draw arrow {7.533000 27.630000 -10.539000} {7.453620 26.979212 -9.783902} red
draw arrow {7.533000 27.630000 -10.539000} {8.528243 27.535336 -10.515961} red
draw arrow {7.533000 27.630000 -10.539000} {7.589488 28.383335 -9.883793} red
draw arrow {7.372000 25.928000 -12.396000} {7.224329 25.050264 -11.940180} blue
draw arrow {7.372000 25.928000 -12.396000} {7.913891 25.470655 -13.101117} blue
draw arrow {7.372000 25.928000 -12.396000} {8.199374 26.070880 -11.852826} blue
