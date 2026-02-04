# Create table to compare CE with SOA
import os
import math
import tabulate
import pandas
import numpy
import matplotlib.pyplot as plt
import re
import collections
import gzip

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
"zeroin.i.3.col":30,
"X-n101-k25.vrp":27591,
"X-n106-k14.vrp":26362,
"X-n110-k13.vrp":14971,
"X-n115-k10.vrp":12747,
"X-n120-k6.vrp":13332,
"X-n125-k30.vrp":55539,
"X-n129-k18.vrp":28940,
"X-n134-k13.vrp":10916,
"X-n139-k10.vrp":13590,
"X-n143-k7.vrp":15700,
"X-n148-k46.vrp":43448,
"X-n153-k22.vrp":21220,
"X-n157-k13.vrp":16876,
"X-n162-k11.vrp":14138,
"X-n167-k10.vrp":20557,
"X-n172-k51.vrp":45607,
"X-n176-k26.vrp":47812,
"X-n181-k23.vrp":25569,
"X-n186-k15.vrp":24145,
"X-n190-k8.vrp":16980,
"X-n195-k51.vrp":44225,
"X-n200-k36.vrp":58578,
"X-n204-k19.vrp":19565,
"X-n209-k16.vrp":30656,
"X-n214-k11.vrp":10856,
"X-n219-k73.vrp":117595,
"X-n223-k34.vrp":40437,
"X-n228-k23.vrp":25742,
"X-n233-k16.vrp":19230,
"X-n237-k14.vrp":27042,
"X-n242-k48.vrp":82751,
"X-n247-k50.vrp":37274,
"X-n251-k28.vrp":38684,
"X-n256-k16.vrp":18839,
"X-n261-k13.vrp":26558,
"X-n266-k58.vrp":75478,
"X-n270-k35.vrp":35291,
"X-n275-k28.vrp":21245,
"X-n280-k17.vrp":33503,
"X-n284-k15.vrp":20215,
"X-n289-k60.vrp":95151,
"X-n294-k50.vrp":47161,
"X-n298-k31.vrp":34231
}

instance_lower_bounds = {
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
"kro124p.1.sop":38762,
"kro124p.2.sop":39841,
"kro124p.3.sop":43904,
"kro124p.4.sop":73021,
"p43.1.sop":28140,
"p43.2.sop":28480,
"p43.3.sop":28835,
"p43.4.sop":83005,
"prob.42.sop":243,
"prob.100.sop":1045,
"rbg048a.sop":351,
"rbg050c.sop":467,
"rbg109a.sop":1038,
"rbg150a.sop":1750,
"rbg174a.sop":2033,
"rbg253a.sop":2950,
"rbg323a.sop":3140,
"rbg341a.sop":2568,
"rbg358a.sop":2545,
"rbg378a.sop":2809,
"ry48p.1.sop":15805,
"ry48p.2.sop":16074,
"ry48p.3.sop":19490,
"ry48p.4.sop":31446,
}

# peel and bound solver
pandb_solver_results_lb = {
"ESC07.sop":2125,
"ESC11.sop":2075,
"ESC12.sop":1675,
"ESC25.sop":1681,
"ESC47.sop":658,
"ESC63.sop":44,
"ESC78.sop":5600,
"br17.10.sop":55,
"br17.12.sop":55,
"ft53.1.sop":4603,
"ft53.2.sop":3555,
"ft53.3.sop":4852,
"ft53.4.sop":7560,
"ft70.1.sop":31122,
"ft70.2.sop":31630,
"ft70.3.sop":32539,
"ft70.4.sop":37984,
"kro124p.1.sop":19224,
"kro124p.2.sop":19299,
"kro124p.3.sop":20145,
"kro124p.4.sop":25002,
"p43.1.sop":27255,
"p43.2.sop":27455,
"p43.3.sop":27780,
"p43.4.sop":28195,
"prob.42.sop":103,
"prob.100.sop":178,
"rbg048a.sop":80,
"rbg050c.sop":175,
"rbg109a.sop":406,
"rbg150a.sop":571,
"rbg174a.sop":646,
"rbg253a.sop":727,
"rbg323a.sop":346,
"rbg341a.sop":340,
"rbg358a.sop":88,
"rbg378a.sop":53,
"ry48p.1.sop":9432,
"ry48p.2.sop":6615,
"ry48p.3.sop":8723,
"ry48p.4.sop":17322
}

pandb_solver_results_time = {
"ESC07.sop":0.04,
"ESC11.sop":0.41,
"ESC12.sop":0.34,
"ESC25.sop":303,
"ESC47.sop":3600,
"ESC63.sop":3600,
"ESC78.sop":3600,
"br17.10.sop":3,
"br17.12.sop":5,
"ft53.1.sop":3600,
"ft53.2.sop":3600,
"ft53.3.sop":3600,
"ft53.4.sop":3600,
"ft70.1.sop":3600,
"ft70.2.sop":3600,
"ft70.3.sop":3600,
"ft70.4.sop":3600,
"kro124p.1.sop":3600,
"kro124p.2.sop":3600,
"kro124p.3.sop":3600,
"kro124p.4.sop":3600,
"p43.1.sop":3600,
"p43.2.sop":3600,
"p43.3.sop":3600,
"p43.4.sop":3600,
"prob.42.sop":3600,
"prob.100.sop":3600,
"rbg048a.sop":3600,
"rbg050c.sop":3600,
"rbg109a.sop":3600,
"rbg150a.sop":3600,
"rbg174a.sop":3600,
"rbg253a.sop":3600,
"rbg323a.sop":3600,
"rbg341a.sop":3600,
"rbg358a.sop":3600,
"rbg378a.sop":3600,
"ry48p.1.sop":3600,
"ry48p.2.sop":3600,
"ry48p.3.sop":3600,
"ry48p.4.sop":3600
}

# Log lines for DDSolver
#STATS - lpIterations[1] lagIterations[4190] sspIterations[4190] numSeparations[4190] numCuts[0] compileTime[0] sspSolveTime[23] lpSolveTime[0] lb[2969.55] ub[83006] numArcs: [165817] numFixed: [0] numHeuristicIPs: [0] numHeuristicLNS: [0] numPrimalLNSRepairs: [0] time: [76]
colelim_pattern = re.compile("STATS - lpIterations\[([0-9]+)\] lagIterations\[([0-9]+)\] sspIterations\[([0-9]+)\] numSeparations\[([0-9]+)\].*compileTime\[([0-9]+)\] sspSolveTime\[([0-9]+)\] lpSolveTime\[([0-9]+)\] lb\[([0-9]+.*)\] ub\[([0-9]+.*)\] numArcs: \[([0-9]+)\] numFixed: \[([0-9]+)\].*time: \[([0-9]+)\]")

logs_dir = "/Users/akarahal/Desktop/dd_vrptw/new_new_logs/"
test_set = ["ALL_sop_LAG_NG2_3600_0"]
#test_set = ["ALL_sop_LAG_NG8"]

instances = []
instance_dir = "/Users/akarahal/Desktop/dd_vrptw/instances/ALL_sop/"

for instance in os.listdir(instance_dir):
  instances.append(instance)

time_results = []
results = []
for test in test_set:
  for instance in instances:
    log_file_name = logs_dir + test + '/' + instance + '.log.gz'
    if not os.path.exists(log_file_name):
      continue

    col_elim = False
    with gzip.open(log_file_name, "rt", encoding='utf-8') as log_file:
      for line in log_file:
        colelim_match = colelim_pattern.match(line)
        if colelim_match:
          col_elim = True
          lpIterations = int(colelim_match.group(1))
          lagIterations = int(colelim_match.group(2))
          sspIterations = colelim_match.group(3)
          numSeparations = int(colelim_match.group(4))
          compileTime = colelim_match.group(5)
          sspSolveTime = colelim_match.group(6)
          lpSolveTime = colelim_match.group(7)
          lb = float(colelim_match.group(8))
          ub = float(colelim_match.group(9))
          numArcs = int(colelim_match.group(10))
          numFixed = int(colelim_match.group(11))
          time = float(colelim_match.group(12))
          time_result = {'instance': instance, 'test': test, 'iterations' : lpIterations, 'lagIterations': lagIterations, 'numSep' : numSeparations, 'lb' : lb, 'ub' : ub, 'numArcs': numArcs, 'numFixed': numFixed, 'time' : time}
          time_results.append(time_result)
   
      if col_elim:
        result = {'instance': instance, 'test': test, 'iterations' : lpIterations, 'lagIterations': lagIterations, 'numSep' : numSeparations, 'lb' : lb, 'numArcs': numArcs, 'numFixed': numFixed, 'time' : time}
        results.append(result)
      else:
        result = {'instance': instance, 'test': test, 'iterations' : numpy.nan, 'lagIterations': numpy.nan, 'numSep' : numpy.nan, 'lb' : numpy.nan, 'numArcs': numpy.nan, 'numFixed': numpy.nan, 'time' : numpy.nan}
        results.append(result)

# cats results
#cats_sol_pattern = re.compile(" \s*[0-9]+ \s*[0-9.]+ \s*[0-9A-Z ]+ \s*[0-9A-Z ]+\s* ([0-9\.]*) \s*([0-9\.]*) \s*BeamSearch.*PE.*")
#cats_sol_pattern = re.compile("\s([0-9]+\.*[0-9]*)\s+([0-9]+\.*[0-9]*)\s+BeamSearch")
cats_sol_pattern = re.compile(" \s*[0-9]+ \s*[0-9.]+ \s*[0-9]+ [A-Z]+ \s*[0-9]+ [A-Z]+ \s* ([0-9\.]*) \s*([0-9\.]*) \s*BeamSearch.*PE.*")
cats_end_dual_pattern = re.compile("^Dual\s*([0-9\.]*)")
cats_end_primal_pattern = re.compile("^Primal\s*([0-9\.]*)")
cats_end_time_pattern = re.compile("^searched\s*([0-9\.]*)")
logs_dir = "/Users/akarahal/Desktop/dd_vrptw/new_new_logs/soa-sop/"
cats_results = {}
for instance in instances:
  log_file_name = logs_dir + '/' + instance + '.log.gz'
  if not os.path.exists(log_file_name):
    continue

  print(log_file_name)
  with gzip.open(log_file_name, "rt", encoding='utf-8') as log_file:
    time = 3600
    lb = -1
    ub = 1e9
    for line in log_file:
      cats_sol_match = cats_sol_pattern.match(line)
      if cats_sol_match:
        lb = cats_sol_match.group(1)
        lb = lb.replace(".","")
        ub = cats_sol_match.group(2)
        ub = ub.replace(".","")
      cats_dual_match = cats_end_dual_pattern.match(line)
      if cats_dual_match:
        lb = cats_dual_match.group(1)
        lb = lb.replace(".","")
      cats_primal_match = cats_end_primal_pattern.match(line)
      if cats_primal_match:
        ub = cats_primal_match.group(1)
        ub = ub.replace(".","")
      cats_time_match = cats_end_time_pattern.match(line)
      if cats_time_match:
        time = cats_time_match.group(1)
        time = int(time.split(".")[0])
    cats_results[instance] = {'lb': lb, 'ub': ub, 'time': time}

print(cats_results)

# give table of results
results_df = pandas.DataFrame(results)
results_df.sort_values(by=['instance','lb'],inplace=True)
results_df.reset_index(drop=True,inplace=True)
table_results_df = results_df[['instance','test','lb','time','iterations','numSep','lagIterations','numArcs']]
print(tabulate.tabulate(table_results_df, headers=table_results_df.columns))

# get sizes for sop instances
pattern = re.compile("DIMENSION: ([0-9]+).*")
instance_sizes = {}
instance_percent_precedence = {}
for instance in instances:
  precedence_count = 0
  instance_file_name = instance_dir + '/' + instance
  instance_file = open(instance_file_name, "r")
  for line in instance_file:
    match = pattern.match(line)
    if match:
      n = int(match.group(1))
      instance_sizes[instance] = n
    split_line = line.split()
    precedence_count = precedence_count + split_line.count('-1')
  instance_percent_precedence[instance] = 100.0 * precedence_count / (n*n)

# create output table for SOA comparison
for i, row in table_results_df.iterrows():
  instance = row['instance']
  instance_name = instance.rsplit('.',1)[0]
  instance_name = instance_name.replace("_","\\_")
  lb_value = row['lb']
  if not math.isnan(lb_value):
    lb_value = round(float(lb_value),1)
  else:
    lb_value = '-'

  numArcs = row['numArcs']
  if not math.isnan(numArcs):
    numArcs = int(numArcs)
  else:
    numArcs = '-'

  numLpIterations = row['iterations']
  if not math.isnan(numLpIterations):
    numLpIterations = int(numLpIterations)
  else:
    numLpIterations = '-'
 
  numLagIterations = row['lagIterations']
  if not math.isnan(numLagIterations):
    numLagIterations = int(numLagIterations)
  else:
    numLagIterations = '-'

  numSeparations = row['numSep']
  if not math.isnan(numSeparations):
    numSeparations = int(numSeparations)
  else:
    numSeparations = '-'

  #gap = round((instance_upper_bounds[instance] - lb_value) * 100.0 / instance_upper_bounds[instance], 1)
  if (lb_value == 0) or (lb_value == '-'):
    lb_value = '-'
    numLpIterations = '-'
    numLagIterations = '-'
    numSeparations = '-'

  time = row['time']
  if not math.isnan(time):
    time = int(time)
    if lb_value != instance_upper_bounds[instance]:
      time = max(3600, time)
  else:
    time = '-'

  percent_precedence = round(instance_percent_precedence[instance],1)

  soa_lb = -1
  soa_ub = 1e9
  soa_time = 3600
  if instance in cats_results:
    soa_lb = cats_results[instance]['lb']
    soa_ub = cats_results[instance]['ub']
    soa_time = cats_results[instance]['time']
  #soa_lb = pandb_solver_results_lb[instance]
  #soa_time = pandb_solver_results_time[instance]

  #print(f"{instance_name} & {instance_upper_bounds[instance]} & & {soa_lb} & {soa_ub} & {soa_time} & & {lb_value} & {numLpIterations} & {numLagIterations} & {numSeparations} & {time} \\\\")
  print(f"{instance_name} & {instance_upper_bounds[instance]} & & {soa_lb} & {soa_time} & & {lb_value} & {numLpIterations} & {numLagIterations} & {numSeparations} & {time} \\\\")
  #print(f"{instance_name} & {instance_sizes[instance]} & {percent_precedence} & & {instance_lower_bounds[instance]} & {instance_upper_bounds[instance]} & & {lb_value} & {numLpIterations} & {numLagIterations} & {numSeparations} & {time} \\\\")
