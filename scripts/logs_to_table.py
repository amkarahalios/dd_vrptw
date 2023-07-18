import os
import math
import tabulate
import pandas
import numpy
import matplotlib.pyplot as plt
import re
import collections

instance_upper_bounds = {
"C101.vrptw":827.3,
"C102.vrptw":827.3,
"C103.vrptw":826.3,
"C104.vrptw":822.9,
"C105.vrptw":827.3,
"C106.vrptw":827.3,
"C107.vrptw":827.3,
"C108.vrptw":827.3,
"C109.vrptw":827.3,
"C201.vrptw":589.1,
"C202.vrptw":589.1,
"C203.vrptw":588.7,
"C204.vrptw":588.1,
"C205.vrptw":586.4,
"C206.vrptw":586,
"C207.vrptw":585.8,
"C208.vrptw":585.8,
"R101.vrptw":1637.7,
"R102.vrptw":1466.6,
"R103.vrptw":1208.7,
"R104.vrptw":971.5,
"R105.vrptw":1355.3,
"R106.vrptw":1234.6,
"R107.vrptw":1064.6,
"R108.vrptw":932.1,
"R109.vrptw":1146.9,
"R110.vrptw":1068,
"R111.vrptw":1048.7,
"R112.vrptw":948.6,
"R201.vrptw":1143.2,
"R202.vrptw":1029.6,
"R203.vrptw":870.8,
"R204.vrptw":731.3,
"R205.vrptw":949.8,
"R206.vrptw":875.9,
"R207.vrptw":794,
"R208.vrptw":701,
"R209.vrptw":854.8,
"R210.vrptw":900.5,
"R211.vrptw":746.7,
"RC101.vrptw":1619.8,
"RC102.vrptw":1457.4,
"RC103.vrptw":1258,
"RC104.vrptw":1132.3,
"RC105.vrptw":1513.7,
"RC106.vrptw":1372.7,
"RC107.vrptw":1207.8,
"RC108.vrptw":1114.2,
"RC201.vrptw":1261.8,
"RC202.vrptw":1092.3,
"RC203.vrptw":923.7,
"RC204.vrptw":783.5,
"RC205.vrptw":1154,
"RC206.vrptw":1051.1,
"RC207.vrptw":962.9,
"RC208.vrptw":776.1,
"C1_2_1.vrptw":2698.6,
"C1_2_2.vrptw":2694.3,
"C1_2_3.vrptw":2675.8,
"C1_2_4.vrptw":2625.6,
"C1_2_5.vrptw":2694.9,
"C1_2_6.vrptw":2694.9,
"C1_2_7.vrptw":2694.9,
"C1_2_8.vrptw":2684,
"C1_2_9.vrptw":2639.6,
"C1_2_10.vrptw":2624.7,
"C2_2_1.vrptw":1922.1,
"C2_2_2.vrptw":1851.4,
"C2_2_3.vrptw":1763.4,
"C2_2_4.vrptw":1695,
"C2_2_5.vrptw":1869.6,
"C2_2_6.vrptw":1844.8,
"C2_2_7.vrptw":1842.2,
"C2_2_8.vrptw":1813.7,
"C2_2_9.vrptw":1815,
"C2_2_10.vrptw":1791.2,
"R1_2_1.vrptw":4667.2,
"R1_2_2.vrptw":3919.9,
"R1_2_3.vrptw":3373.9,
"R1_2_4.vrptw":3047.6,
"R1_2_5.vrptw":4053.2,
"R1_2_6.vrptw":3559.1,
"R1_2_7.vrptw":3141.9,
"R1_2_8.vrptw":2938.4,
"R1_2_9.vrptw":3734.7,
"R1_2_10.vrptw":3293.1,
"R2_2_1.vrptw":3468,
"R2_2_2.vrptw":3008.2,
"R2_2_3.vrptw":2537.5,
"R2_2_4.vrptw":1928.5,
"R2_2_5.vrptw":3061.1,
"R2_2_6.vrptw":2675.4,
"R2_2_7.vrptw":2304.7,
"R2_2_8.vrptw":1842.4,
"R2_2_9.vrptw":2843.3,
"R2_2_10.vrptw":2549.4,
"RC1_2_1.vrptw":3516.9,
"RC1_2_2.vrptw":3221.6,
"RC1_2_3.vrptw":3001.4,
"RC1_2_4.vrptw":2845.2,
"RC1_2_5.vrptw":3325.6,
"RC1_2_6.vrptw":3300.7,
"RC1_2_7.vrptw":3177.8,
"RC1_2_8.vrptw":3060,
"RC1_2_9.vrptw":3073.3,
"RC1_2_10.vrptw":2990.5,
"RC2_2_1.vrptw":2797.4,
"RC2_2_2.vrptw":2481.6,
"RC2_2_3.vrptw":2227.7,
"RC2_2_4.vrptw":1854.8,
"RC2_2_5.vrptw":2491.4,
"RC2_2_6.vrptw":2495.1,
"RC2_2_7.vrptw":2287.7,
"RC2_2_8.vrptw":2151.2,
"RC2_2_9.vrptw":2086.6,
"RC2_2_10.vrptw":1989.2,
"C1_4_1.vrptw":7138.8,
"C1_4_2.vrptw":7113.3,
"C1_4_3.vrptw":6929.9,
"C1_4_4.vrptw":6777.7,
"C1_4_5.vrptw":7138.8,
"C1_4_6.vrptw":7140.1,
"C1_4_7.vrptw":7136.2,
"C1_4_8.vrptw":7083,
"C1_4_9.vrptw":6927.8,
"C1_4_10.vrptw":6825.4,
"C2_4_1.vrptw":4100.3,
"C2_4_2.vrptw":3914.1,
"C2_4_3.vrptw":3755.2,
"C2_4_4.vrptw":3523.7,
"C2_4_5.vrptw":3923.2,
"C2_4_6.vrptw":3860.1,
"C2_4_7.vrptw":3870.9,
"C2_4_8.vrptw":3773.7,
"C2_4_9.vrptw":3842.1,
"C2_4_10.vrptw":3665.1,
"R1_4_1.vrptw":10305.8,
"R1_4_2.vrptw":8873.3,
"R1_4_3.vrptw":7784.3,
"R1_4_4.vrptw":7266.2,
"R1_4_5.vrptw":9184.6,
"R1_4_6.vrptw":8340.4,
"R1_4_7.vrptw":7599.8,
"R1_4_8.vrptw":7240.5,
"R1_4_9.vrptw":8677.5,
"R1_4_10.vrptw":8077.8,
"R2_4_1.vrptw":7520.7,
"R2_4_2.vrptw":6482.8,
"R2_4_3.vrptw":5372.9,
"R2_4_4.vrptw":4211.2,
"R2_4_5.vrptw":6567.9,
"R2_4_6.vrptw":5813.5,
"R2_4_7.vrptw":4893.5,
"R2_4_8.vrptw":4000.1,
"R2_4_9.vrptw":6067.8,
"R2_4_10.vrptw":5645.9,
"RC1_4_1.vrptw":8522.9,
"RC1_4_2.vrptw":7878.2,
"RC1_4_3.vrptw":7516.9,
"RC1_4_4.vrptw":7292.9,
"RC1_4_5.vrptw":8152.3,
"RC1_4_6.vrptw":8148,
"RC1_4_7.vrptw":7932.5,
"RC1_4_8.vrptw":7757.2,
"RC1_4_9.vrptw":7717.7,
"RC1_4_10.vrptw":7581.2,
"RC2_4_1.vrptw":6147.3,
"RC2_4_2.vrptw":5407.5,
"RC2_4_3.vrptw":4573,
"RC2_4_4.vrptw":3597.9,
"RC2_4_5.vrptw":5392.3,
"RC2_4_6.vrptw":5324.6,
"RC2_4_7.vrptw":4987.8,
"RC2_4_8.vrptw":4693.3,
"RC2_4_9.vrptw":4510.4,
"RC2_4_10.vrptw":4252.3,
"C1_6_1.vrptw":14076.6,
"C1_6_2.vrptw":13948.3,
"C1_6_3.vrptw":13757,
"C1_6_4.vrptw":13538.6,
"C1_6_5.vrptw":14066.8,
"C1_6_6.vrptw":14070.9,
"C1_6_7.vrptw":14066.8,
"C1_6_8.vrptw":13991.2,
"C1_6_9.vrptw":13664.5,
"C1_6_10.vrptw":13617.5,
"C2_6_1.vrptw":7752.2,
"C2_6_2.vrptw":7471.5,
"C2_6_3.vrptw":7215,
"C2_6_4.vrptw":6877,
"C2_6_5.vrptw":7553.8,
"C2_6_6.vrptw":7449.8,
"C2_6_7.vrptw":7491.3,
"C2_6_8.vrptw":7303.7,
"C2_6_9.vrptw":7303.2,
"C2_6_10.vrptw":7123.9,
"R1_6_1.vrptw":21274.2,
"R1_6_2.vrptw":18558.7,
"R1_6_3.vrptw":16874.9,
"R1_6_4.vrptw":15721.4,
"R1_6_5.vrptw":19294.9,
"R1_6_6.vrptw":17763.7,
"R1_6_7.vrptw":16496.2,
"R1_6_8.vrptw":15584.3,
"R1_6_9.vrptw":18474.1,
"R1_6_10.vrptw":17583.7,
"R2_6_1.vrptw":15145.3,
"R2_6_2.vrptw":12976.3,
"R2_6_3.vrptw":10455.3,
"R2_6_4.vrptw":7915.1,
"R2_6_5.vrptw":13790.2,
"R2_6_6.vrptw":11847.8,
"R2_6_7.vrptw":9777.9,
"R2_6_8.vrptw":7512.3,
"R2_6_9.vrptw":12736.8,
"R2_6_10.vrptw":11837,
"RC1_6_1.vrptw":16960.1,
"RC1_6_2.vrptw":15890.6,
"RC1_6_3.vrptw":15181.3,
"RC1_6_4.vrptw":14753.2,
"RC1_6_5.vrptw":16536.3,
"RC1_6_6.vrptw":16473.3,
"RC1_6_7.vrptw":16055.3,
"RC1_6_8.vrptw":15891.8,
"RC1_6_9.vrptw":15803.5,
"RC1_6_10.vrptw":15651.3,
"RC2_6_1.vrptw":11966.1,
"RC2_6_2.vrptw":10336.9,
"RC2_6_3.vrptw":8894.9,
"RC2_6_4.vrptw":6967.5,
"RC2_6_5.vrptw":11080.7,
"RC2_6_6.vrptw":10830.5,
"RC2_6_7.vrptw":10289.4,
"RC2_6_8.vrptw":9779,
"RC2_6_9.vrptw":9436,
"RC2_6_10.vrptw":8974.7,
"C1_8_1.vrptw":25156.9,
"C1_8_2.vrptw":24974.1,
"C1_8_3.vrptw":24156.1,
"C1_8_4.vrptw":23797.3,
"C1_8_5.vrptw":25138.6,
"C1_8_6.vrptw":25133.3,
"C1_8_7.vrptw":25127.3,
"C1_8_8.vrptw":24809.7,
"C1_8_9.vrptw":24200.4,
"C1_8_10.vrptw":24026.7,
"C2_8_1.vrptw":11631.9,
"C2_8_2.vrptw":11394.5,
"C2_8_3.vrptw":11138.1,
"C2_8_4.vrptw":10650,
"C2_8_5.vrptw":11395.6,
"C2_8_6.vrptw":11316.3,
"C2_8_7.vrptw":11332.9,
"C2_8_8.vrptw":11133.9,
"C2_8_9.vrptw":11140.4,
"C2_8_10.vrptw":10946,
"R1_8_1.vrptw":36345,
"R1_8_2.vrptw":32277.6,
"R1_8_3.vrptw":29304.5,
"R1_8_4.vrptw":27734.7,
"R1_8_5.vrptw":33494.2,
"R1_8_6.vrptw":30872.4,
"R1_8_7.vrptw":28789,
"R1_8_8.vrptw":27609.4,
"R1_8_9.vrptw":32257.3,
"R1_8_10.vrptw":30918.4,
"R2_8_1.vrptw":24969.8,
"R2_8_2.vrptw":21312.2,
"R2_8_3.vrptw":17234.8,
"R2_8_4.vrptw":13160.8,
"R2_8_5.vrptw":22801.6,
"R2_8_6.vrptw":19740.5,
"R2_8_7.vrptw":16357.5,
"R2_8_8.vrptw":12611.7,
"R2_8_9.vrptw":21282.7,
"R2_8_10.vrptw":19984.8,
"RC1_8_1.vrptw":29978.9,
"RC1_8_2.vrptw":28290.1,
"RC1_8_3.vrptw":27447.7,
"RC1_8_4.vrptw":26557.2,
"RC1_8_5.vrptw":29219.9,
"RC1_8_6.vrptw":29194.2,
"RC1_8_7.vrptw":28788.6,
"RC1_8_8.vrptw":28418.1,
"RC1_8_9.vrptw":28347.1,
"RC1_8_10.vrptw":28168.5,
"RC2_8_1.vrptw":19201.3,
"RC2_8_2.vrptw":16709.5,
"RC2_8_3.vrptw":14013.6,
"RC2_8_4.vrptw":10969.4,
"RC2_8_5.vrptw":17466.1,
"RC2_8_6.vrptw":17195.1,
"RC2_8_7.vrptw":16362.2,
"RC2_8_8.vrptw":15528.8,
"RC2_8_9.vrptw":15183,
"RC2_8_10.vrptw":14370.9,
"C1_10_1.vrptw":42444.8,
"C1_10_2.vrptw":41337.8,
"C1_10_3.vrptw":40064.4,
"C1_10_4.vrptw":39434.1,
"C1_10_5.vrptw":42434.8,
"C1_10_6.vrptw":42437,
"C1_10_7.vrptw":42420.4,
"C1_10_8.vrptw":41652.1,
"C1_10_9.vrptw":40288.4,
"C1_10_10.vrptw":39816.8,
"C2_10_1.vrptw":16841.1,
"C2_10_2.vrptw":16462.6,
"C2_10_3.vrptw":16036.5,
"C2_10_4.vrptw":15459.5,
"C2_10_5.vrptw":16521.3,
"C2_10_6.vrptw":16290.7,
"C2_10_7.vrptw":16378.4,
"C2_10_8.vrptw":16029.1,
"C2_10_9.vrptw":16075.4,
"C2_10_10.vrptw":15728.6,
"R1_10_1.vrptw":53046.5,
"R1_10_2.vrptw":48263.1,
"R1_10_3.vrptw":44677.1,
"R1_10_4.vrptw":42440.7,
"R1_10_5.vrptw":50406.7,
"R1_10_6.vrptw":46930.3,
"R1_10_7.vrptw":43997.4,
"R1_10_8.vrptw":42279.3,
"R1_10_9.vrptw":49162.8,
"R1_10_10.vrptw":47364.6,
"R2_10_1.vrptw":36881,
"R2_10_2.vrptw":31241.9,
"R2_10_3.vrptw":24399,
"R2_10_4.vrptw":17811.5,
"R2_10_5.vrptw":34132.8,
"R2_10_6.vrptw":29124.7,
"R2_10_7.vrptw":23102.2,
"R2_10_8.vrptw":17403.8,
"R2_10_9.vrptw":31990.6,
"R2_10_10.vrptw":29840.5,
"RC1_10_1.vrptw":45790.8,
"RC1_10_2.vrptw":43678.3,
"RC1_10_3.vrptw":42122,
"RC1_10_4.vrptw":41357.4,
"RC1_10_5.vrptw":45028.1,
"RC1_10_6.vrptw":44903.6,
"RC1_10_7.vrptw":44417.1,
"RC1_10_8.vrptw":43916.5,
"RC1_10_9.vrptw":43858.1,
"RC1_10_10.vrptw":43533.7,
"RC2_10_1.vrptw":28122.6,
"RC2_10_2.vrptw":24248.6,
"RC2_10_3.vrptw":19618.1,
"RC2_10_4.vrptw":15657,
"RC2_10_5.vrptw":25797.5,
"RC2_10_6.vrptw":25782.5,
"RC2_10_7.vrptw":24395.8,
"RC2_10_8.vrptw":23280.2,
"RC2_10_9.vrptw":22731.6,
"RC2_10_10.vrptw":21736.1,
"ESC07.sop":2125,
"ESC11.sop":2075,
"ESC12.sop":1675,
"ESC25.sop":1681,
"ESC47.sop":1288,
"ESC63.sop":62,
"ESC78.sop":18230,
"br17.10.sop":55,
"br17.12.sop":55,
"ft53.1.sop":7531,
"ft53.2.sop":8026,
"ft53.3.sop":10262,
"ft53.4.sop":14425,
"ft70.1.sop":39313,
"ft70.2.sop":40419,
"ft70.3.sop":42535,
"ft70.4.sop":53530,
"kro124p.1.sop":39420,
"kro124p.2.sop":41336,
"kro124p.3.sop":49499,
"kro124p.4.sop":76103,
"p43.1.sop":28140,
"p43.2.sop":28480,
"p43.3.sop":28835,
"p43.4.sop":83005,
"prob.42.sop":243,
"prob.100.sop":1163,
"rbg048a.sop":351,
"rbg050c.sop":467,
"rbg109a.sop":1038,
"rbg150a.sop":1750,
"rbg174a.sop":2033,
"rbg253a.sop":2950,
"rbg323a.sop":3140,
"rbg341a.sop":2568,
"rbg358a.sop":2545,
"rbg378a.sop":2816,
"ry48p.1.sop":15805,
"ry48p.2.sop":16666,
"ry48p.3.sop":19894,
"ry48p.4.sop":31446,
"rc201.0.tsptw":628.62,
"rc201.1.tsptw":654.7,
"rc201.2.tsptw":707.65,
"rc201.3.tsptw":422.54,
"rc202.0.tsptw":496.22,
"rc202.1.tsptw":426.53,
"rc202.2.tsptw":611.77,
"rc202.3.tsptw":627.85,
"rc203.0.tsptw":727.45,
"rc203.1.tsptw":726.99,
"rc203.2.tsptw":617.46,
"rc204.0.tsptw":541.45,
"rc204.1.tsptw":485.37,
"rc204.2.tsptw":778.4,
"rc205.0.tsptw":511.65,
"rc205.1.tsptw":491.22,
"rc205.2.tsptw":714.69,
"rc205.3.tsptw":601.24,
"rc206.0.tsptw":835.23,
"rc206.1.tsptw":664.73,
"rc206.2.tsptw":655.37,
"rc207.0.tsptw":806.69,
"rc207.1.tsptw":726.36,
"rc207.2.tsptw":546.41,
"rc208.0.tsptw":820.56,
"rc208.1.tsptw":509.04,
"rc208.2.tsptw":503.92,
"rc_201.1.tsptw":444.54,
"rc_201.2.tsptw":711.54,
"rc_201.3.tsptw":790.61,
"rc_201.4.tsptw":793.64,
"rc_202.1.tsptw":771.78,
"rc_202.2.tsptw":304.14,
"rc_202.3.tsptw":837.72,
"rc_202.4.tsptw":793.03,
"rc_203.1.tsptw":453.48,
"rc_203.2.tsptw":784.16,
"rc_203.3.tsptw":817.53,
"rc_203.4.tsptw":314.29,
"rc_204.1.tsptw":878.64,
"rc_204.2.tsptw":662.16,
"rc_204.3.tsptw":455.03,
"rc_205.1.tsptw":343.21,
"rc_205.2.tsptw":755.93,
"rc_205.3.tsptw":825.06,
"rc_205.4.tsptw":760.47,
"rc_206.1.tsptw":117.85,
"rc_206.2.tsptw":828.06,
"rc_206.3.tsptw":574.42,
"rc_206.4.tsptw":831.67,
"rc_207.1.tsptw":732.68,
"rc_207.2.tsptw":701.25,
"rc_207.3.tsptw":682.4,
"rc_207.4.tsptw":119.64,
"rc_208.1.tsptw":789.25,
"rc_208.2.tsptw":533.78,
"rc_208.3.tsptw":634.44,
"1-FullIns_3.col":4,
"1-FullIns_4.col":5,
"1-FullIns_5.col":6,
"1-Insertions_4.col":5,
"1-Insertions_5.col":6,
"1-Insertions_6.col":7,
"2-FullIns_3.col":5,
"2-FullIns_4.col":6,
"2-FullIns_5.col":7,
"2-Insertions_3.col":4,
"2-Insertions_4.col":5,
"2-Insertions_5.col":6,
"3-FullIns_3.col":6,
"3-FullIns_4.col":7,
"3-FullIns_5.col":8,
"3-Insertions_3.col":4,
"3-Insertions_4.col":5,
"3-Insertions_5.col":6,
"4-FullIns_3.col":7,
"4-FullIns_4.col":8,
"4-FullIns_5.col":9,
"4-Insertions_3.col":4,
"4-Insertions_4.col":5,
"5-FullIns_3.col":8,
"5-FullIns_4.col":9,
"abb313GPIA.col":9,
"anna.col":11,
"ash331GPIA.col":4,
"ash608GPIA.col":4,
"ash958GPIA.col":4,
"C2000.5.col":145,
"C2000.9.col":400,
"C4000.5.col":259,
"david.col":11,
"DSJC1000.1.col":20,
"DSJC1000.5.col":82,
"DSJC1000.9.col":222,
"DSJC125.1.col":5,
"DSJC125.5.col":17,
"DSJC125.9.col":44,
"DSJC250.1.col":8,
"DSJC250.5.col":28,
"DSJC250.9.col":72,
"DSJC500.1.col":12,
"DSJC500.5.col":47,
"DSJC500.9.col":126,
"DSJR500.1.col":12,
"DSJR500.1c.col":85,
"DSJR500.5.col":122,
"flat1000_50_0.col":50,
"flat1000_60_0.col":60,
"flat1000_76_0.col":81,
"flat300_20_0.col":20,
"flat300_26_0.col":26,
"flat300_28_0.col":28,
"fpsol2.i.1.col":65,
"fpsol2.i.2.col":30,
"fpsol2.i.3.col":30,
"games120.col":9,
"homer.col":13,
"huck.col":11,
"inithx.i.1.col":54,
"inithx.i.2.col":31,
"inithx.i.3.col":31,
"jean.col":10,
"latin_square_10.col":98,
"le450_15a.col":15,
"le450_15b.col":15,
"le450_15c.col":15,
"le450_15d.col":15,
"le450_25a.col":25,
"le450_25b.col":25,
"le450_25c.col":25,
"le450_25d.col":25,
"le450_5a.col":5,
"le450_5b.col":5,
"le450_5c.col":5,
"le450_5d.col":5,
"miles1000.col":42,
"miles1500.col":73,
"miles250.col":8,
"miles500.col":20,
"miles750.col":31,
"mug100_1.col":4,
"mug100_25.col":4,
"mug88_1.col":4,
"mug88_25.col":4,
"mulsol.i.1.col":49,
"mulsol.i.2.col":31,
"mulsol.i.3.col":31,
"mulsol.i.4.col":31,
"mulsol.i.5.col":31,
"myciel3.col":4,
"myciel4.col":5,
"myciel5.col":6,
"myciel6.col":7,
"myciel7.col":8,
"qg.order100.col":100,
"qg.order30.col":30,
"qg.order40.col":40,
"qg.order60.col":60,
"queen10_10.col":11,
"queen11_11.col":11,
"queen12_12.col":12,
"queen13_13.col":13,
"queen14_14.col":14,
"queen15_15.col":15,
"queen16_16.col":16,
"queen5_5.col":5,
"queen6_6.col":7,
"queen7_7.col":7,
"queen8_12.col":12,
"queen8_8.col":9,
"queen9_9.col":10,
"r1000.1.col":20,
"r1000.1c.col":97,
"r1000.5.col":234,
"r125.1.col":5,
"r125.1c.col":46,
"r125.5.col":36,
"r250.1.col":8,
"r250.1c.col":64,
"r250.5.col":65,
"school1.col":14,
"school1_nsh.col":14,
"wap01a.col":41,
"wap02a.col":41,
"wap03a.col":44,
"wap04a.col":42,
"wap05a.col":50,
"wap06a.col":40,
"wap07a.col":41,
"wap08a.col":41,
"will199GPIA.col":7,
"zeroin.i.1.col":49,
"zeroin.i.2.col":30,
"zeroin.i.3.col":30
}

soa_bounds = {
}

# Log lines for DDSolver
#STATS - lpIterations[1] lagIterations[32] sspIterations[186] numSeparations[0] compileTime[240] sspSolveTime[349] lpSolveTime[0] lb[794.025] ub[1e+10] size: [87254] time: [590]

# Log lines for ColGenSolver
#STATS - iterations[402] ddTime[50] pricingTime[112] masterTime[97] lb[1581.24] time: [260]

colelim_pattern = re.compile("STATS - lpIterations\[([0-9]+)\] lagIterations\[([0-9]+)\] sspIterations\[([0-9]+)\] numSeparations\[([0-9]+)\].*compileTime\[([0-9]+)\] sspSolveTime\[([0-9]+)\] lpSolveTime\[([0-9]+)\] lb\[([0-9]+.*)\] ub\[([0-9]+.*)\] size: \[([0-9]+)\] time: \[([0-9]+)\]")
colgen_pattern = re.compile("STATS - iterations\[([0-9]+)\] ddTime\[([0-9]+)\] pricingTime\[([0-9]+)\] masterTime\[([0-9]+)\] lb\[([0-9]+.*)\] time: \[([0-9]+)\]")
colgen_pattern1 = re.compile("STATS - iterations\[([0-9]+)\] ddTime\[([0-9]+)\] pricingTime\[([0-9]+)\] masterTime\[([0-9]+)\] lb\[([0-9]+.*)\] sol\[([0-9]+.*)\] time: \[([0-9]+)\]")
colgen_complete_pattern = re.compile("no negative reduced cost paths")
colgen_pattern_size = re.compile("DD size: ([0-9]+)")
held_initial_pattern = re.compile("Finished initial bounds: LB ([0-9]+) and UB ([0-9]+) in.*([0-9]+)\.0000.*seconds.")
held_update_time_pattern = re.compile("Compute_coloring took ([0-9]+)..* seconds.*")
held_improved_lb_pattern = re.compile("Lower bound improved: LB ([0-9]+) and UB ([0-9]+).*")
held_improved_ub_pattern = re.compile("Upper bound improved: LB ([0-9]+) and UB ([0-9]+).*")
bdd_pattern = re.compile(".*Time elapsed: ([0-9]+).*BDD.*Lower bound: ([0-9]+).*Upper bound: ([0-9]+).*")

root_node_pattern = re.compile(".*Time elapsed: ([0-9]+).*BDD.*Lower bound: ([0-9]+).*Upper bound: ([0-9]+).*")
time_pattern = re.compile("Current Bounds - LB.*UB.* Time: ([0-9]+).*")
lb_pattern = re.compile("Current Bounds - LB: ([0-9]+) UB:.*")
ub_pattern = re.compile("Current Bounds - LB:.*UB: ([0-9]+).*")
finish_pattern = re.compile("Done solving.*time:\[([0-9]+).*\] LB:\[([0-9]+)\] UB:\[([0-9]+)\]")

#logs_dir = "/Users/akarahal/Desktop/dd_vrptw/logs/"
logs_dir = "/Users/akarahal/Desktop/dd_graph_color/logs/"
#test_set = ["col_elim_hg_lag_ng_4_20_N_3600","col_elim_hg-close_lag_ng_4_20_N-7200"]
#test_set = ["col_elim_hg-close_lag_ng_4_20_N-7200"]
test_set = ["col_elim_hg_lag_ng_2_50_N_3600",
            "col_elim_hg_lag_ng_4_50_N_3600",
            "col_elim_hg_lag_ng_6_50_N_3600",
            "col_elim_hg_lp_ng_2_50_N_3600",
            "col_elim_hg_lp_ng_4_50_N_3600",
            "col_elim_hg_lp_ng_6_50_N_3600",
            "col_elim_sop_lag_ng_2_50_N_3600",
            "col_elim_sop_lag_ng_5_50_N_3600",
            "col_elim_sop_lag_ng_8_50_N_3600",
            "col_elim_sop_lp_ng_2_50_N_3600",
            "col_elim_sop_lp_ng_5_50_N_3600",
            "col_elim_sop_lp_ng_8_50_N_3600",
            "col_elim_spb_lag_ng_2_50_N_3600",
            "col_elim_spb_lag_ng_5_50_N_3600",
            "col_elim_spb_lag_ng_8_50_N_3600",
            "col_elim_spb_lp_ng_2_50_N_3600",
            "col_elim_spb_lp_ng_5_50_N_3600",
            "col_elim_spb_lp_ng_8_50_N_3600",
            "col_elim_lp_0_0_BANDB_3600",
            "col_elim_lp_5_0_BANDB_3600",
            "col_elim_lp_0_0_MIP_3600",
            "col_elim_lp_5_0_MIP_3600",
            "col_elim_lp_5_3_MIP_3600"]

instances = []
#instance_dir = "/Users/akarahal/Desktop/dd_vrptw/instances/Vrp-Set-HG/"
#instance_dir = "/Users/akarahal/Desktop/dd_vrptw/instances/ALL_sop/"
#instance_dir = "/Users/akarahal/Desktop/dd_vrptw/instances/SolomonPotvinBengio/"
instance_dir = "/Users/akarahal/Desktop/dd_graph_color/instances/"

for instance in os.listdir(instance_dir):
  if "COG" in instance:
    continue
  instances.append(instance)

time_results = []
results = []
for test in test_set:
  for instance in instances:
    log_file_name = logs_dir + test + '/' + instance + '.log'
    if not os.path.exists(log_file_name):
      continue

    col_elim = False
    col_gen = False
    col_gen_complete = False
    held = False
    log_file = open(log_file_name, "r")
    for line in log_file:
      colelim_match = colelim_pattern.match(line)
      if colelim_match:
        col_elim = True
        lpIterations = colelim_match.group(1)
        lagIterations = colelim_match.group(2)
        sspIterations = colelim_match.group(3)
        numSeparations = colelim_match.group(4)
        compileTime = colelim_match.group(5)
        sspSolveTime = colelim_match.group(6)
        lpSolveTime = colelim_match.group(7)
        lb = math.ceil(float(colelim_match.group(8)))
        ub = float(colelim_match.group(9))
        size = colelim_match.group(10)
        time = colelim_match.group(11)
        time_result = {'instance': instance, 'test': test, 'iterations' : lpIterations, 'numSep' : numSeparations, 'lb' : lb, 'ub' : ub, 'size': size, 'time' : time}
        time_results.append(time_result)
 
      colgen_match = colgen_pattern.match(line)
      colgen_match1 = colgen_pattern1.match(line)
      if colgen_match and not colgen_match1:
        col_gen = True
        iterations = colgen_match.group(1)
        ddTime = colgen_match.group(2)
        pricingTime = colgen_match.group(3)
        masterTime = colgen_match.group(4)
        lb = math.ceil(float(colgen_match.group(5)))
        time = colgen_match.group(6)
        # remove DD time?
        time = float(time) - float(ddTime)
 
      if colgen_match1:
        col_gen = True
        iterations = colgen_match1.group(1)
        ddTime = colgen_match1.group(2)
        pricingTime = colgen_match1.group(3)
        masterTime = colgen_match1.group(4)
        lb = math.ceil(float(colgen_match1.group(5)))
        sol = float(colgen_match1.group(6))
        time = colgen_match1.group(7)
        # remove DD time?
        time = float(time) - float(ddTime)
        time_result = {'instance': instance, 'test': test, 'iterations' : iterations, 'lb' : lb, 'size' : size, 'time' : time}
        time_results.append(time_result)

      colgen_size_match = colgen_pattern_size.match(line)
      if colgen_size_match:
        size = colgen_size_match.group(1)

      colgen_complete_match = colgen_complete_pattern.match(line)
      if colgen_complete_match:
        col_gen_complete = True
 
      held_initial_match = held_initial_pattern.match(line)
      if held_initial_match:
        held = True
        lb = float(held_initial_match.group(1))
        ub = float(held_initial_match.group(2))
        time = float(held_initial_match.group(3))
        time_result = {'instance': instance, 'test': test, 'lb' : lb, 'ub' : ub, 'time' : time}
        time_results.append(time_result)
 
      held_update_time_match = held_update_time_pattern.match(line)
      if held_update_time_match:
        time = float(held_update_time_match.group(1))

      held_improved_lb_match = held_improved_lb_pattern.match(line)
      if held_improved_lb_match:
        lb = float(held_improved_lb_match.group(1))
        ub = float(held_improved_lb_match.group(2))
        time_result = {'instance': instance, 'test': test, 'lb' : lb, 'ub' : ub, 'time' : time}
        time_results.append(time_result)
 
      held_improved_ub_match = held_improved_ub_pattern.match(line)
      if held_improved_ub_match:
        lb = float(held_improved_ub_match.group(1))
        ub = float(held_improved_ub_match.group(2))
        time_result = {'instance': instance, 'test': test, 'lb' : lb, 'ub' : ub, 'time' : time}
        time_results.append(time_result)

      bdd_match = bdd_pattern.match(line)
      if bdd_match:
        time = float(bdd_match.group(1))
        lb = float(bdd_match.group(2))
        ub = float(bdd_match.group(3))
        time_result = {'instance': instance, 'test': test, 'lb' : lb, 'ub' : ub, 'time' : time}
        time_results.append(time_result)

    if col_gen and col_gen_complete:
      result = {'instance': instance, 'test': test, 'iterations' : iterations, 'lb' : lb, 'size': size, 'time' : time}
      results.append(result)
      time_results.append(result)
    elif col_elim:
      result = {'instance': instance, 'test': test, 'iterations' : lpIterations, 'numSep' : numSeparations, 'lb' : lb, 'size': size, 'time' : time}
      results.append(result)
    elif held:
      result = {'instance': instance, 'test': test, 'lb' : lb, 'time' : time}
      results.append(result)
    else:
      result = {'instance': instance, 'test': test, 'iterations' : numpy.nan, 'numSep' : numpy.nan, 'lb' : numpy.nan, 'size': numpy.nan, 'time' : numpy.nan}
      results.append(result)

# setup results df
results_df = pandas.DataFrame(results)
results_df.sort_values(by=['instance','lb'],inplace=True)
results_df.reset_index(drop=True,inplace=True)
time_results_df = pandas.DataFrame(time_results)
time_results_df.sort_values(by=['instance','lb'],inplace=True)
time_results_df.reset_index(drop=True,inplace=True)

# print table
print_table = True
if print_table:
  table_results_df = results_df[['instance','test','lb','time']]
  print(tabulate.tabulate(table_results_df, headers=table_results_df.columns))

# experiment 1 - Exact DD Size vs. DD Elim Size
experiment1 = False
if experiment1:
  x = []
  y = []
  label = []
  for instance in set(results_df['instance']):
    instance_results = results_df[results_df['instance'] == instance]
    instance_tests = set(instance_results['test'])
    #if (("test_A_colelim_lag_q_2_2_1800" in instance_tests) and ("test_A_colelim_lp_q_1_2_1800" in instance_tests) and ("test_A_colgen_dd_q_2_1800" in instance_tests)):
    #  colelim_q_2_2 = instance_results[instance_results['test'] == "test_A_colelim_lag_q_2_2_1800"].iloc[0]
    #  colelim_q_1_2 = instance_results[instance_results['test'] == "test_A_colelim_lp_q_1_2_1800"].iloc[0]
    #  colgen_q_2_2 = instance_results[instance_results['test'] == "test_A_colgen_dd_q_2_1800"].iloc[0]
    #if (("A_colelim_lp_q_1_2_1800" in instance_tests) and ("A_colelim_lp_q_2_2_1800" in instance_tests) and ("A_colgen_dd_q_2_1800" in instance_tests)):
    if (("A_colelim_lp_q_1_2_1800" in instance_tests) and ("A_colgen_dd_q_2_1800" in instance_tests)):
      colelim_q_1_2 = instance_results[instance_results['test'] == "A_colelim_lp_q_1_2_1800"].iloc[0]
      colgen_q_2_2 = instance_results[instance_results['test'] == "A_colgen_dd_q_2_1800"].iloc[0]
      if (colgen_q_2_2['lb'] == colelim_q_1_2['lb']):
        exact_size = int(colgen_q_2_2['size'])
        elim_size = int(colelim_q_1_2['size'])
        x.append(exact_size)
        y.append(elim_size)
        label.append(instance)
        plt.plot([exact_size],[elim_size],marker='o',label=instance)
    if (("B_colelim_lp_q_1_2_1800" in instance_tests) and ("B_colgen_dd_q_2_1800" in instance_tests)):
      colelim_q_1_2 = instance_results[instance_results['test'] == "B_colelim_lp_q_1_2_1800"].iloc[0]
      colgen_q_2_2 = instance_results[instance_results['test'] == "B_colgen_dd_q_2_1800"].iloc[0]
      if (colgen_q_2_2['lb'] == colelim_q_1_2['lb']):
        exact_size = int(colgen_q_2_2['size'])
        elim_size = int(colelim_q_1_2['size'])
        x.append(exact_size)
        y.append(elim_size)
        label.append(instance)
        plt.plot([exact_size],[elim_size],marker='o',label=instance)
  plt.xlabel('exact size')
  plt.ylabel('elim size at opt')
  plt.xscale('log')
  plt.yscale('log')
  plt.xlim([10**2,10**5])
  plt.ylim([10**2,10**5])
  #plt.legend()
  plt.title('q=1->2 vs. q=2')
  plt.show()

# experiment 2 - plot the lb versus time for each method
experiment2 = False
if experiment2:
  # used these for paper exp2
  # using these to test X instances
  #tests_to_compare = ["col_elim_hg_lag_ng_4_20_N_3600","col_elim_hg-close_lag_ng_4_20_N-7200"]
  tests_to_compare = ["col_elim_sop_lag_ng_2_50_N_3600",
                      "col_elim_sop_lag_ng_5_50_N_3600",
                      "col_elim_sop_lag_ng_8_50_N_3600",
                      "col_elim_sop_lp_ng_2_50_N_3600",
                      "col_elim_sop_lp_ng_5_50_N_3600",
                      "col_elim_sop_lp_ng_8_50_N_3600"]

  #instances_to_consider = ["X-n327-k20.vrp","X-n344-k43.vrp","X-n359-k29.vrp","X-n367-k17.vrp","X-n480-k70.vrp","X-n502-k39.vrp","X-n561-k42.vrp","X-n573-k30.vrp","X-n801-k40.vrp","X-n957-k87.vrp"]
  instances_to_consider = instances
  #instances_to_consider = instances['X']
  for instance in instances_to_consider:
    #time = list(range(0,3600,100))
    time = list(range(0,7200,100))
    lbs = []
    labels = []
    instance_results = time_results_df[time_results_df['instance'] == instance]
    if instance_results.empty:
      continue
    if not instance_results.empty:
      for test in tests_to_compare:
        instance_test_results = instance_results[instance_results['test'] == test]
        if not instance_test_results.empty:
          instance_test_lb = []
          for t in time:
            instance_test_results_time = instance_test_results[instance_test_results['time'].astype(int) <= t]
            if (instance_test_results_time.empty):
              lb = 0
            else:
              lb = max(instance_test_results_time['lb'])
            instance_test_lb.append(lb)
          lbs.append(instance_test_lb)
          labels.append(test)

      min_val = 10000000
      max_val = 0
      for (label, lb) in zip(labels, lbs):
        if len([l for l in lb if l >0]) == 0:
          continue
        to_plot = zip(time,lb)
        to_plot = [data for data in to_plot if data[1] > 0]
        times_to_plot = [d[0] for d in to_plot]
        lbs_to_plot = [d[1] for d in to_plot]
        plt.plot(times_to_plot,lbs_to_plot,label=label)
        min_value = min([l for l in lb if l > 0])
        max_value = max(lb)
        min_val = min(min_val, min_value)
        max_val = max(max_val, max_value)
      plt.title(instance)
      plt.xlabel('time (s)')
      plt.ylabel('lb')
      if instance in instance_upper_bounds:
        plt.axhline(y=instance_upper_bounds[instance],linewidth=2,label="opt",color="magenta")
        plt.ylim(min_val-10,instance_upper_bounds[instance]+10)
      else:
        plt.ylim(min_val,max_val)
      plt.legend()
      plt.show()

# experiment 3 - lp runtimes
experiment3 = False
if experiment3:
  tests_to_include = ["test_A_colelim_lp_ng_2_10_1800", "test_A_colelim_lp_ng_5_10_1800","test_M_colgen_dd_q_1_1800"]
  instances_to_consider = ["A-n32-k5.vrp", "A-n33-k5.vrp","A-n36-k5.vrp","A-n37-k6.vrp","A-n38-k5.vrp","A-n80-k10.vrp","M-n101-k10.vrp","M-n121-k7.vrp","M-n151-k12.vrp"]
  #instances_to_consider = ["A-n80-k10.vrp","M-n101-k10.vrp","M-n121-k7.vrp","M-n151-k12.vrp"]
  ddSizes = []
  lpSolveTimes = []
  labels = []
  for instance in instances_to_consider:
    instance_results = time_results_df[time_results_df['instance'] == instance]
    for test in tests_to_include:
      instance_test_results = instance_results[instance_results['test'] == test]
      if not instance_test_results.empty:
        instance_test_lb = []
        curr_time = 0
        instance_test_results.dropna(inplace=True)
        for index, row in instance_test_results.iterrows():
          solve_time = int(row['time']) - curr_time
          curr_time = int(row['time'])
          dd_size = int(row['size'])
          ddSizes.append(dd_size)
          lpSolveTimes.append(solve_time)
          if 'M' in instance:
            labels.append('blue')
          elif 'A' in instance:
            labels.append('orange')

  plt.scatter(ddSizes, lpSolveTimes, color=labels)
  plt.title('LP time vs. DD size')
  plt.xlabel('DD size')
  plt.ylabel('time (s)')
  plt.show()


'''                   
                    'CELAGmuSSP' : ['test_A_colelim_lag_q_1_2_N_3600',
                                    'test_B_colelim_lag_q_1_2_N_3600',
                                    'test_M_colelim_lag_q_1_2_N_3600'],
'''
# experiment Baseline
# We want this to be a performance plot
# First showing how many instances are solved within d% at times 0-3600
# Then showing how many instances are solved within x% of the optimal value
experiment4 = False
if experiment4:
  methodsToTests = {'CE_MIP' : ['all_col_elim_lp_ub_3600'],
                    'CE_BANDB' : ['all_col_elim_lp_bandb_3600'],
                    'Held' : ['held_3600']}
  # d% gap
  # List of times
  # List of percent gaps
  # x_axis names?
  d = 5
  times_list = list(range(0,3700,100))
  gap_list = list(range(d+1,21,1))
  line_styles = ['solid', 'dotted', 'dashed', 'dashdot']
  line_colors = ['b', 'r', 'g', 'k']

  # For each method build:
  # - num solved within d% at each time
  # - num solved within each percent gaps
  method_num = 0
  for method in methodsToTests:
    count = 0
    method_time_list = [0] * len(times_list)
    method_gap_list = [0] * len(gap_list)
    for test in methodsToTests[method]:
      for instance in instances:
        count = count + 1
        optimal = instance_upper_bounds[instance]
        test_instance_results = time_results_df[(time_results_df['instance'] == instance) & (time_results_df['test'] == test)]
        if not test_instance_results.empty:
          print(method)
          print(test)
          # time info
          for t_iter in range(len(times_list)):
            t = times_list[t_iter]
            time_result = test_instance_results[test_instance_results['time'].astype(int) <= t]
            if time_result.empty:
              lb = 0
            else:
              lb = max(time_result['lb'])
            if lb >= (optimal * (1-(d*1.0/100))):
              method_time_list[t_iter] = method_time_list[t_iter] + 1
          # gap info
          best_lb = max(test_instance_results['lb'])
          for gap_iter in range(len(gap_list)):
            gap = gap_list[gap_iter]
            if best_lb >= (optimal * (1- (0.01*gap))):
              method_gap_list[gap_iter] = method_gap_list[gap_iter] + 1

    # plot each of these time then reverse gaps?
    x_plot = times_list.copy()
    for gap_iter in range(len(gap_list)):
        x_plot.append(3600 + (3600 / len(gap_list) * (gap_iter+1)))
    y_plot = method_time_list + method_gap_list
    plt.plot(x_plot,y_plot,label=method,linestyle=line_styles[method_num],color=line_colors[method_num])
    method_num = method_num + 1
    print("method count")
    print(count)

  # plot
  plt.xlabel('time (s) | optimality gap (%)')
  plt.ylabel('# instances')

  # update x axis
  x_ticks = [100,1000,1800,3600]
  x_labels = [100,1000,1800,3600]
  #x_ticks_raw = times_list.copy()
  #x_labels_raw = times_list.copy()

  for gap_iter in range(len(gap_list)):
    if gap_iter % 2 == 0 and gap_iter > 0:
      x_ticks.append(3600 + (3600 / len(gap_list) * (gap_iter+1)))
      x_labels.append(str(gap_list[gap_iter]) + '%')

  plt.xticks(x_ticks, x_labels)
  plt.legend()
  plt.axvline(x=3600,color='k')
  plt.show()

# Get the average optimality gap for each class of instances
experiment5 = False
if experiment5:
  instance_type_tests = {'A' : ['test_A_colelim_lp_ng_2_20_Y_3600', 'test_A_colelim_lag_ng_2_20_Y_3600'],
                         'B' : ['test_B_colelim_lp_ng_2_20_Y_3600', 'test_B_colelim_lag_ng_2_20_Y_3600'],
                         'E' : ['test_E_colelim_lp_ng_2_20_Y_3600', 'test_E_colelim_lag_ng_2_20_Y_3600'],
                         'F' : ['test_F_colelim_lp_ng_2_20_Y_3600', 'test_F_colelim_lag_ng_2_20_Y_3600'],
                         'M' : ['test_M_colelim_lp_ng_2_20_Y_3600', 'test_M_colelim_lag_ng_2_20_Y_3600'],
                         'P' : ['test_P_colelim_lp_ng_2_20_Y_3600', 'test_P_colelim_lag_ng_2_20_Y_3600'],
                         'X' : ['test_X_colelim_lp_ng_2_20_Y_7200', 'test_X_colelim_lag_ng_2_20_Y_7200']}

  best_gaps = {}
  instance_set_list = ['A','B','E','F','M','P','X']
  for instance_type in instance_set_list:
    instance_type_gaps = []
    for instance in instances[instance_type]:
      # Remove instances not in Pecin
      q2_non_compile = ['A-n32','A-n33','A-n34','A-n36','B-n31','B-n34','B-n35']
      shouldSkip = False
      for q2_no in q2_non_compile:
        if q2_no in instance:
          shouldSkip = True
          break
      if shouldSkip:
        continue

      optimal = instance_optimal_values[instance]
      instance_best_gap = 10000000000
      for test in instance_type_tests[instance_type]:
        test_instance_results = time_results_df[(time_results_df['instance'] == instance) & (time_results_df['test'] == test)]
        if not test_instance_results.empty:
          lb = max(test_instance_results['lb'])
        else:
          lb = 0

        gap = (optimal - lb) * 100.0 / optimal
        instance_best_gap = min(instance_best_gap, gap)
      if instance_best_gap != 100.0:
        instance_type_gaps.append(instance_best_gap)
      else:
        print(instance)
    print(instance_type_gaps)
    instance_type_key = instance_type
    if instance_type == 'M':
      instance_type_key = 'E'
    best_gaps[instance_type_key] = 1.0 * sum(instance_type_gaps) / len(instance_type_gaps)
  print(best_gaps)

# experiment Baseline
# We want this to be a performance plot
# Don't use best ub like experiment4 though
# First showing how many instances are solved within d% at times 0-3600
# Then showing how many instances are solved within x% of the optimal value
experiment6 = True
if experiment6:
  #methodsToTests = {'CE_MIP' : ['all_col_elim_lp_ub_3600'],
  #                  'CE_BANDB' : ['all_col_elim_lp_bandb_ub_3600'],
  #                  'CE_BDD' : ['all_col_elim_bdd_3600'],
  #                  'Held' : ['held_3600']}
  #methodsToTests = {'LP_NG2' : ['col_elim_sop_lp_ng_2_50_N_3600'],
  #                  'LP_NG5' : ['col_elim_sop_lp_ng_5_50_N_3600'],
  #                  'LP_NG8' : ['col_elim_sop_lp_ng_8_50_N_3600'],
  #                  'LAG_NG2' : ['col_elim_sop_lag_ng_2_50_N_3600'],
  #                  'LAG_NG5' : ['col_elim_sop_lag_ng_5_50_N_3600'],
  #                  'LAG_NG8' : ['col_elim_sop_lag_ng_8_50_N_3600']}
  #methodsToTests = {'LP_NG2' : ['col_elim_spb_lp_ng_2_50_N_3600'],
  #                  'LP_NG5' : ['col_elim_spb_lp_ng_5_50_N_3600'],
  #                  'LP_NG8' : ['col_elim_spb_lp_ng_8_50_N_3600'],
  #                  'LAG_NG2' : ['col_elim_spb_lag_ng_2_50_N_3600'],
  #                  'LAG_NG5' : ['col_elim_spb_lag_ng_5_50_N_3600'],
  #                  'LAG_NG8' : ['col_elim_spb_lag_ng_8_50_N_3600']}
  #methodsToTests = {'LP_NG2' : ['col_elim_hg_lp_ng_2_50_N_3600'],
  #                  'LP_NG4' : ['col_elim_hg_lp_ng_4_50_N_3600'],
  #                  'LP_NG6' : ['col_elim_hg_lp_ng_6_50_N_3600'],
  #                  'LAG_NG2' : ['col_elim_hg_lag_ng_2_50_N_3600'],
  #                  'LAG_NG4' : ['col_elim_hg_lag_ng_4_50_N_3600'],
  #                  'LAG_NG6' : ['col_elim_hg_lag_ng_6_50_N_3600']}
  methodsToTests = {'MIP_0_0' : ['col_elim_lp_0_0_MIP_3600'],
                    'MIP_5_0' : ['col_elim_lp_5_0_MIP_3600'],
                    'MIP_5_3' : ['col_elim_lp_5_3_MIP_3600'],
                    'BANDB_0_0' : ['col_elim_lp_0_0_BANDB_3600'],
                    'BANDB_5_0' : ['col_elim_lp_5_0_BANDB_3600']}

  # List of times
  # List of percent gaps
  # x_axis names?
  d = 0
  times_list = list(range(0,3700,100))
  gap_list = list(range(d,21,1))
  line_styles = ['solid', 'dotted', 'dashed', 'dashdot', (0, (5, 1)), (5, (10, 3))]
  line_colors = ['b', 'r', 'g', 'k', 'saddlebrown', 'm']

  # For each method build:
  # - num solved within d% at each time
  # - num solved within each percent gaps
  method_num = 0
  for method in methodsToTests:
    count = 0
    method_time_list = [0] * len(times_list)
    method_gap_list = [0] * len(gap_list)
    for test in methodsToTests[method]:
      for instance in instances:
        count = count + 1
        optimal = instance_upper_bounds[instance]
        test_instance_results = time_results_df[(time_results_df['instance'] == instance) & (time_results_df['test'] == test)]
        if not test_instance_results.empty:
          print(method)
          print(test)
          # time info
          for t_iter in range(len(times_list)):
            t = times_list[t_iter]
            time_result = test_instance_results[test_instance_results['time'].astype(int) <= t]
            if time_result.empty:
              lb = 0
              ub = 1000000
            else:
              lb = max(time_result['lb'])
              ub = min(time_result['ub'])
            if lb >= ub * (1-(d*1.0/100)):
              method_time_list[t_iter] = method_time_list[t_iter] + 1
          # gap info
          best_lb = max(test_instance_results['lb'])
          best_ub = min(test_instance_results['ub'])
          for gap_iter in range(len(gap_list)):
            gap = gap_list[gap_iter]
            if best_lb >= (best_ub * (1- (0.01*gap))):
              method_gap_list[gap_iter] = method_gap_list[gap_iter] + 1

    # plot each of these time then reverse gaps?
    x_plot = times_list.copy()
    for gap_iter in range(len(gap_list)):
        x_plot.append(3600 + (3600 / len(gap_list) * (gap_iter+1)))
    y_plot = method_time_list + method_gap_list
    plt.plot(x_plot,y_plot,label=method,linestyle=line_styles[method_num],color=line_colors[method_num])
    method_num = method_num + 1
    print("method count")
    print(count)

  # plot
  plt.xlabel('time (s) | optimality gap (%)')
  plt.ylabel('# instances')

  # update x axis
  x_ticks = [100,1000,1800,3600]
  x_labels = [100,1000,1800,3600]
  #x_ticks_raw = times_list.copy()
  #x_labels_raw = times_list.copy()

  for gap_iter in range(len(gap_list)):
    if gap_iter % 2 == 0 and gap_iter > 0:
      x_ticks.append(3600 + (3600 / len(gap_list) * (gap_iter+1)))
      x_labels.append(str(gap_list[gap_iter]) + '%')

  plt.xticks(x_ticks, x_labels)
  plt.legend()
  plt.axvline(x=3600,color='k')
  plt.show()
