mol new 2jqx.THR.roundtrip.xyz
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
topo addbond 0 14
topo addbond 1 0
topo addbond 2 1
topo addbond 2 15
topo addbond 3 2
topo addbond 4 1
topo addbond 5 4
topo addbond 6 4
topo addbond 7 0
topo addbond 8 1
topo addbond 9 4
topo addbond 10 5
topo addbond 11 6
topo addbond 12 6
topo addbond 13 6
draw arrow {14.849000 27.195000 -12.966000} {15.210067 28.036832 -12.564815} red
draw arrow {14.849000 27.195000 -12.966000} {14.684855 27.675863 -13.827294} red
draw arrow {14.849000 27.195000 -12.966000} {13.931020 27.440132 -12.654194} red
draw arrow {16.546000 26.094000 -11.628000} {16.087275 26.083414 -12.516515} blue
draw arrow {16.546000 26.094000 -11.628000} {17.088916 25.299082 -11.898827} blue
draw arrow {16.546000 26.094000 -11.628000} {15.842570 25.487376 -11.257604} blue
