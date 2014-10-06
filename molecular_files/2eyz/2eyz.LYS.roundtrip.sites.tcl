mol new 2eyz.LYS.roundtrip.xyz
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
color Name L blue3
color Name M blue
color Name P blue
color Name Q green2
color Name V green
color Name W green
color Name X black
proc vmd_draw_arrow {mol start end color} {
    # an arrow is made of a cylinder and a cone
    set middle [vecadd $start [vecscale 0.800000 [vecsub $end $start]]]
    graphics $mol color $color
    graphics $mol cylinder $start $middle radius 0.100000
    graphics $mol cone $middle $end radius 0.200000
}
topo addbond 0 22
topo addbond 1 0
topo addbond 2 1
topo addbond 2 23
topo addbond 3 2
topo addbond 4 1
topo addbond 5 4
topo addbond 6 5
topo addbond 7 6
topo addbond 8 7
topo addbond 9 0
topo addbond 10 1
topo addbond 11 4
topo addbond 12 4
topo addbond 13 5
topo addbond 14 5
topo addbond 15 6
topo addbond 16 6
topo addbond 17 7
topo addbond 18 7
topo addbond 19 8
topo addbond 20 8
topo addbond 21 8
draw arrow {462.638000 -4.512000 3.173000} {462.755032 -5.373304 3.667428} red
draw arrow {462.638000 -4.512000 3.173000} {461.870689 -4.274350 3.768615} red
draw arrow {462.638000 -4.512000 3.173000} {462.007494 -4.961086 2.539925} red
draw arrow {463.893000 -2.393000 3.926000} {463.932541 -3.269503 3.446230} blue
draw arrow {463.893000 -2.393000 3.926000} {464.752591 -2.118360 3.495099} blue
draw arrow {463.893000 -2.393000 3.926000} {464.402450 -2.788368 4.690294} blue
draw arrow {464.985000 -0.969000 2.166000} {464.045992 -0.722494 1.926210} green
draw arrow {464.985000 -0.969000 2.166000} {464.675757 -1.269226 3.068349} green
draw arrow {464.985000 -0.969000 2.166000} {465.135444 -0.047534 2.524145} green
