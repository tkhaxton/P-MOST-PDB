mol new 2eyz.Nterm.MET.roundtrip.xyz
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
topo addbond 1 0
topo addbond 2 1
topo addbond 2 19
topo addbond 3 2
topo addbond 4 1
topo addbond 5 4
topo addbond 6 5
topo addbond 7 6
topo addbond 8 1
topo addbond 9 4
topo addbond 10 4
topo addbond 11 5
topo addbond 12 5
topo addbond 13 7
topo addbond 14 7
topo addbond 15 7
topo addbond 16 0
topo addbond 17 0
topo addbond 18 0
draw arrow {474.007000 27.891000 -10.910000} {473.647525 28.780503 -11.192069} red
draw arrow {474.007000 27.891000 -10.910000} {473.144407 27.689559 -10.445937} red
draw arrow {474.007000 27.891000 -10.910000} {474.362965 28.301130 -10.070308} red
draw arrow {476.879000 25.968000 -11.637000} {476.725455 25.146676 -12.186410} blue
draw arrow {476.879000 25.968000 -11.637000} {476.269457 26.484342 -12.238539} blue
draw arrow {476.879000 25.968000 -11.637000} {477.656741 26.210526 -12.216913} blue
