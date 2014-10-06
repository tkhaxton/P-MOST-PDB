mol new 2jqx.CYS.roundtrip.xyz
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
draw arrow {-26.302000 5.559000 17.949000} {-25.755758 5.096795 18.647560} red
draw arrow {-26.302000 5.559000 17.949000} {-27.022403 5.725265 18.622332} red
draw arrow {-26.302000 5.559000 17.949000} {-26.729363 4.687953 17.706847} red
draw arrow {-24.977000 6.483000 15.656000} {-25.287244 6.627596 16.595596} blue
draw arrow {-24.977000 6.483000 15.656000} {-24.953994 5.496066 15.815477} blue
draw arrow {-24.977000 6.483000 15.656000} {-24.026621 6.554093 15.958864} blue
