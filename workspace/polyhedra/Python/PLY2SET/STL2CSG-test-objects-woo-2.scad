module s0()
{
 translate([-0.0, 0.0, 0.0])
 rotate(a = -90.0, v = [0.0, 1.0, -0.0])
 translate([-103.07764064044152, -103.07764064044152, -206.15528128088303])
 cube([206.15528128088303, 206.15528128088303, 206.15528128088303]);
}

module s2()
{
 translate([0.0, 0.0, 10.0])
 translate([-103.07764064044152, -103.07764064044152, -206.15528128088303])
 cube([206.15528128088303, 206.15528128088303, 206.15528128088303]);
}

module s3()
{
 translate([0.0, 0.0, 10.0])
 translate([-103.07764064044152, -103.07764064044152, 0.0])
 cube([206.15528128088303, 206.15528128088303, 206.15528128088303]);
}

module s4()
{
 translate([0.0, -0.0, 0.0])
 rotate(a = -90.0, v = [-1.0, 0.0, 0.0])
 translate([-103.07764064044152, -103.07764064044152, -206.15528128088303])
 cube([206.15528128088303, 206.15528128088303, 206.15528128088303]);
}

module s6()
{
 translate([0.0, 0.0, -0.0])
 translate([-103.07764064044152, -103.07764064044152, 0.0])
 cube([206.15528128088303, 206.15528128088303, 206.15528128088303]);
}

module s8()
{
 translate([10.0, 0.0, 0.0])
 rotate(a = -90.0, v = [0.0, -1.0, 0.0])
 translate([-103.07764064044152, -103.07764064044152, -206.15528128088303])
 cube([206.15528128088303, 206.15528128088303, 206.15528128088303]);
}

module s10()
{
 translate([0.0, 10.0, 0.0])
 rotate(a = -90.0, v = [1.0, 0.0, 0.0])
 translate([-103.07764064044152, -103.07764064044152, -206.15528128088303])
 cube([206.15528128088303, 206.15528128088303, 206.15528128088303]);
}

module s12()
{
 translate([0.0, 8.0, 0.0])
 rotate(a = -90.0, v = [1.0, 0.0, 0.0])
 translate([-103.07764064044152, -103.07764064044152, -206.15528128088303])
 cube([206.15528128088303, 206.15528128088303, 206.15528128088303]);
}

module s14()
{
 translate([2.0, -0.0, -0.0])
 rotate(a = -90.0, v = [0.0, 1.0, -0.0])
 translate([-103.07764064044152, -103.07764064044152, -206.15528128088303])
 cube([206.15528128088303, 206.15528128088303, 206.15528128088303]);
}

module s16()
{
 translate([0.0, 0.0, 15.0])
 translate([-103.07764064044152, -103.07764064044152, -206.15528128088303])
 cube([206.15528128088303, 206.15528128088303, 206.15528128088303]);
}

module s18()
{
 translate([-0.0, 2.0, -0.0])
 rotate(a = -90.0, v = [-1.0, 0.0, 0.0])
 translate([-103.07764064044152, -103.07764064044152, -206.15528128088303])
 cube([206.15528128088303, 206.15528128088303, 206.15528128088303]);
}

module s20()
{
 translate([8.0, 0.0, 0.0])
 rotate(a = -90.0, v = [0.0, -1.0, 0.0])
 translate([-103.07764064044152, -103.07764064044152, -206.15528128088303])
 cube([206.15528128088303, 206.15528128088303, 206.15528128088303]);
}

module s23()
{
 translate([-0.0, -3.4482763, 1.3793105])
 rotate(a = -68.19859976473145, v = [-0.9284767508506775, 0.0, 0.0])
 translate([-103.07764064044152, -103.07764064044152, -206.15528128088303])
 cube([206.15528128088303, 206.15528128088303, 206.15528128088303]);
}

module s24()
{
 translate([-3.4482763, 0.0, 1.3793106])
 rotate(a = -68.198590569306, v = [0.0, 0.9284766912460327, -0.0])
 translate([-103.07764064044152, -103.07764064044152, -206.15528128088303])
 cube([206.15528128088303, 206.15528128088303, 206.15528128088303]);
}

module s26()
{
 translate([0.0, 12.068966, 4.8275867])
 rotate(a = -68.198590569306, v = [0.9284766912460327, 0.0, 0.0])
 translate([-103.07764064044152, -103.07764064044152, -206.15528128088303])
 cube([206.15528128088303, 206.15528128088303, 206.15528128088303]);
}

module s28()
{
 translate([12.068968, 0.0, 4.827587])
 rotate(a = -68.19859976473145, v = [0.0, -0.9284767508506775, 0.0])
 translate([-103.07764064044152, -103.07764064044152, -206.15528128088303])
 cube([206.15528128088303, 206.15528128088303, 206.15528128088303]);
}

union(){intersection(){s0();intersection(){s2();intersection(){s4();intersection(){s6();intersection(){s8();intersection(){s10();intersection(){s16();intersection(){s23();intersection(){s24();intersection(){s26();s28();};};};};};};};};};};intersection(){s3();intersection(){s12();intersection(){s14();intersection(){s16();intersection(){s18();s20();};};};};};}
