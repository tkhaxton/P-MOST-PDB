mol new 2eyz.MET.roundtrip.xyz
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
color Name I blue2
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
topo addbond 0 17
topo addbond 1 0
topo addbond 2 1
topo addbond 2 18
topo addbond 3 2
topo addbond 4 1
topo addbond 5 4
topo addbond 6 5
topo addbond 7 6
topo addbond 8 0
topo addbond 9 1
topo addbond 10 4
topo addbond 11 4
topo addbond 12 5
topo addbond 13 5
topo addbond 14 7
topo addbond 15 7
topo addbond 16 7
draw arrow {457.966000 2.793000 -10.953000} {457.803667 1.856966 -10.640770} red
draw arrow {457.966000 2.793000 -10.953000} {457.108438 2.770311 -11.466880} red
draw arrow {457.966000 2.793000 -10.953000} {458.454093 2.441824 -11.752025} red
draw arrow {456.518000 5.238000 -8.089000} {455.835948 4.969871 -7.408624} blue
draw arrow {456.518000 5.238000 -8.089000} {456.665627 4.276311 -8.320000} blue
draw arrow {456.518000 5.238000 -8.089000} {457.234248 5.180888 -7.393495} blue
