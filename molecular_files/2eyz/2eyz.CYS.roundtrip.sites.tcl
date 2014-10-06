mol new 2eyz.CYS.roundtrip.xyz
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
topo addbond 0 11
topo addbond 1 0
topo addbond 2 1
topo addbond 2 12
topo addbond 3 2
topo addbond 4 1
topo addbond 5 4
topo addbond 6 0
topo addbond 7 1
topo addbond 8 4
topo addbond 9 4
topo addbond 10 5
draw arrow {486.687000 -3.386000 8.179000} {487.406117 -4.055160 7.991661} red
draw arrow {486.687000 -3.386000 8.179000} {485.992788 -4.089718 8.027832} red
draw arrow {486.687000 -3.386000 8.179000} {486.656322 -3.147240 7.208406} red
draw arrow {487.928000 -1.314000 9.606000} {487.347376 -2.127980 9.588305} blue
draw arrow {487.928000 -1.314000 9.606000} {487.131976 -0.751015 9.828248} blue
draw arrow {487.928000 -1.314000 9.606000} {487.757057 -1.170872 8.631170} blue
