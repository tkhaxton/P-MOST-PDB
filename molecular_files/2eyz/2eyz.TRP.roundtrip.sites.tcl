mol new 2eyz.TRP.roundtrip.xyz
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
topo addbond 0 24
topo addbond 1 0
topo addbond 2 1
topo addbond 2 25
topo addbond 3 2
topo addbond 4 1
topo addbond 5 4
topo addbond 6 5
topo addbond 7 5
topo addbond 8 6
topo addbond 9 7
topo addbond 9 8
topo addbond 10 7
topo addbond 11 9
topo addbond 12 10
topo addbond 13 11
topo addbond 13 12
topo addbond 14 0
topo addbond 15 1
topo addbond 16 4
topo addbond 17 4
topo addbond 18 6
topo addbond 19 8
topo addbond 20 10
topo addbond 21 11
topo addbond 22 12
topo addbond 23 13
draw arrow {460.686000 10.689000 -1.241000} {460.862985 11.398252 -0.558623} red
draw arrow {460.686000 10.689000 -1.241000} {461.603394 10.821225 -1.616373} red
draw arrow {460.686000 10.689000 -1.241000} {460.329538 11.381444 -1.868262} red
draw arrow {461.455000 8.868000 0.308000} {461.682777 9.156564 1.237972} blue
draw arrow {461.455000 8.868000 0.308000} {462.226945 8.232363 0.316162} blue
draw arrow {461.455000 8.868000 0.308000} {462.048480 9.584028 -0.059539} blue
