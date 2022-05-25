\ Model JSuperPoly
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  0 K_16 + 0 flag_16 + 0 K_17 + 0 flag_17 + 0 K_18 + 0 flag_18 + 0 K_19
   + 0 flag_19 + 0 K_20 + 0 flag_20 + 0 K_21 + 0 flag_21 + 0 K_22
   + 0 flag_22 + 0 K_23 + 0 flag_23 + 0 K_24 + 0 flag_24 + 0 K_25
   + 0 flag_25 + 0 K_26 + 0 flag_26 + 0 K_27 + 0 flag_27 + 0 K_28
   + 0 flag_28 + 0 K_29 + 0 flag_29 + 0 K_30 + 0 flag_30 + 0 K_31
   + 0 flag_31 + 0 K_32 + 0 flag_32 + 0 K_33 + 0 flag_33 + 0 K_34
   + 0 flag_34 + 0 K_35 + 0 flag_35 + 0 K_36 + 0 flag_36 + 0 K_37
   + 0 flag_37 + 0 K_38 + 0 flag_38 + 0 K_39 + 0 flag_39 + 0 K_40
   + 0 flag_40 + 0 K_41 + 0 flag_41 + 0 K_42 + 0 flag_42 + 0 K_43
   + 0 flag_43 + 0 K_44 + 0 flag_44 + 0 K_45 + 0 flag_45 + 0 K_46
   + 0 flag_46 + 0 K_47 + 0 flag_47 + 0 K_48 + 0 flag_48 + 0 K_49
   + 0 flag_49 + 0 K_50 + 0 flag_50 + 0 K_51 + 0 flag_51 + 0 K_52
   + 0 flag_52 + 0 K_53 + 0 flag_53 + 0 K_54 + 0 flag_54 + 0 K_55
   + 0 flag_55 + 0 K_56 + 0 flag_56 + 0 K_57 + 0 flag_57 + 0 K_58
   + 0 flag_58 + 0 K_59 + 0 flag_59 + 0 K_60 + 0 flag_60 + 0 K_61
   + 0 flag_61 + 0 K_62 + 0 flag_62 + 0 K_63 + 0 flag_63
Subject To
 R0: - RK^0_0 + RK^0_0_copy1 <= 0
 R1: - RK^0_0 + RK^0_0_copy2 <= 0
 R2: RK^0_0 - RK^0_0_copy1 - RK^0_0_copy2 <= 0
 R3: - L^0_0 + L^0_0_copy1 <= 0
 R4: - L^0_0 + L^0_0_copy2 <= 0
 R5: L^0_0 - L^0_0_copy1 - L^0_0_copy2 <= 0
 R6: a^1_0 - RK^0_0_copy2 = 0
 R7: a^1_0 - L^0_0_copy2 = 0
 R8: - b^1_0 + R^0_0_copy2 = 0
 R9: - b^1_0 + c^1_0 >= 0
 R10: - a^1_0 - b^1_0 + c^1_0 = 0
 R11: - c^1_1 + c^1_1_0_copy1 <= 0
 R12: - c^1_1 + c^1_1_0_copy2 <= 0
 R13: c^1_1 - c^1_1_0_copy1 - c^1_1_0_copy2 <= 0
 R14: - c^1_2 + c^1_2_0_copy1 <= 0
 R15: - c^1_2 + c^1_2_0_copy2 <= 0
 R16: c^1_2 - c^1_2_0_copy1 - c^1_2_0_copy2 <= 0
 R17: d^1_0 - c^1_1_0_copy2 = 0
 R18: d^1_0 - c^1_2_0_copy2 = 0
 R19: e^1_0 - d^1_0 >= 0
 R20: - c^1_0 + e^1_0 - d^1_0 = 0
 R21: - RK^0_1 + RK^0_1_copy1 <= 0
 R22: - RK^0_1 + RK^0_1_copy2 <= 0
 R23: RK^0_1 - RK^0_1_copy1 - RK^0_1_copy2 <= 0
 R24: - L^0_1 + L^0_1_copy1 <= 0
 R25: - L^0_1 + L^0_1_copy2 <= 0
 R26: L^0_1 - L^0_1_copy1 - L^0_1_copy2 <= 0
 R27: a^1_1 - RK^0_1_copy2 = 0
 R28: a^1_1 - L^0_1_copy2 = 0
 R29: - b^1_1 + R^0_1_copy2 = 0
 R30: - b^1_1 + c^1_1_0_copy1 >= 0
 R31: - a^1_1 - b^1_1 + c^1_1_0_copy1 = 0
 R32: - c^1_2_0_copy1 + c^1_1_1_copy1 <= 0
 R33: - c^1_2_0_copy1 + c^1_1_1_copy2 <= 0
 R34: c^1_2_0_copy1 - c^1_1_1_copy1 - c^1_1_1_copy2 <= 0
 R35: - c^1_3 + c^1_2_1_copy1 <= 0
 R36: - c^1_3 + c^1_2_1_copy2 <= 0
 R37: c^1_3 - c^1_2_1_copy1 - c^1_2_1_copy2 <= 0
 R38: d^1_1 - c^1_1_1_copy2 = 0
 R39: d^1_1 - c^1_2_1_copy2 = 0
 R40: e^1_1 - d^1_1 >= 0
 R41: e^1_1 - d^1_1 - c^1_1_0_copy1 = 0
 R42: - RK^0_2 + RK^0_2_copy1 <= 0
 R43: - RK^0_2 + RK^0_2_copy2 <= 0
 R44: RK^0_2 - RK^0_2_copy1 - RK^0_2_copy2 <= 0
 R45: - L^0_2 + L^0_2_copy1 <= 0
 R46: - L^0_2 + L^0_2_copy2 <= 0
 R47: L^0_2 - L^0_2_copy1 - L^0_2_copy2 <= 0
 R48: a^1_2 - RK^0_2_copy2 = 0
 R49: a^1_2 - L^0_2_copy2 = 0
 R50: - b^1_2 + R^0_2_copy2 = 0
 R51: - b^1_2 + c^1_1_1_copy1 >= 0
 R52: - a^1_2 - b^1_2 + c^1_1_1_copy1 = 0
 R53: - c^1_2_1_copy1 + c^1_1_2_copy1 <= 0
 R54: - c^1_2_1_copy1 + c^1_1_2_copy2 <= 0
 R55: c^1_2_1_copy1 - c^1_1_2_copy1 - c^1_1_2_copy2 <= 0
 R56: - c^1_4 + c^1_2_2_copy1 <= 0
 R57: - c^1_4 + c^1_2_2_copy2 <= 0
 R58: c^1_4 - c^1_2_2_copy1 - c^1_2_2_copy2 <= 0
 R59: d^1_2 - c^1_1_2_copy2 = 0
 R60: d^1_2 - c^1_2_2_copy2 = 0
 R61: e^1_2 - d^1_2 >= 0
 R62: e^1_2 - d^1_2 - c^1_1_1_copy1 = 0
 R63: - RK^0_3 + RK^0_3_copy1 <= 0
 R64: - RK^0_3 + RK^0_3_copy2 <= 0
 R65: RK^0_3 - RK^0_3_copy1 - RK^0_3_copy2 <= 0
 R66: - L^0_3 + L^0_3_copy1 <= 0
 R67: - L^0_3 + L^0_3_copy2 <= 0
 R68: L^0_3 - L^0_3_copy1 - L^0_3_copy2 <= 0
 R69: a^1_3 - RK^0_3_copy2 = 0
 R70: a^1_3 - L^0_3_copy2 = 0
 R71: - b^1_3 + R^0_3_copy2 = 0
 R72: - b^1_3 + c^1_1_2_copy1 >= 0
 R73: - a^1_3 - b^1_3 + c^1_1_2_copy1 = 0
 R74: - c^1_2_2_copy1 + c^1_1_3_copy1 <= 0
 R75: - c^1_2_2_copy1 + c^1_1_3_copy2 <= 0
 R76: c^1_2_2_copy1 - c^1_1_3_copy1 - c^1_1_3_copy2 <= 0
 R77: - c^1_5 + c^1_2_3_copy1 <= 0
 R78: - c^1_5 + c^1_2_3_copy2 <= 0
 R79: c^1_5 - c^1_2_3_copy1 - c^1_2_3_copy2 <= 0
 R80: d^1_3 - c^1_1_3_copy2 = 0
 R81: d^1_3 - c^1_2_3_copy2 = 0
 R82: e^1_3 - d^1_3 >= 0
 R83: e^1_3 - d^1_3 - c^1_1_2_copy1 = 0
 R84: - RK^0_4 + RK^0_4_copy1 <= 0
 R85: - RK^0_4 + RK^0_4_copy2 <= 0
 R86: RK^0_4 - RK^0_4_copy1 - RK^0_4_copy2 <= 0
 R87: - L^0_4 + L^0_4_copy1 <= 0
 R88: - L^0_4 + L^0_4_copy2 <= 0
 R89: L^0_4 - L^0_4_copy1 - L^0_4_copy2 <= 0
 R90: a^1_4 - RK^0_4_copy2 = 0
 R91: a^1_4 - L^0_4_copy2 = 0
 R92: - b^1_4 + R^0_4_copy2 = 0
 R93: - b^1_4 + c^1_1_3_copy1 >= 0
 R94: - a^1_4 - b^1_4 + c^1_1_3_copy1 = 0
 R95: - c^1_2_3_copy1 + c^1_1_4_copy1 <= 0
 R96: - c^1_2_3_copy1 + c^1_1_4_copy2 <= 0
 R97: c^1_2_3_copy1 - c^1_1_4_copy1 - c^1_1_4_copy2 <= 0
 R98: - c^1_6 + c^1_2_4_copy1 <= 0
 R99: - c^1_6 + c^1_2_4_copy2 <= 0
 R100: c^1_6 - c^1_2_4_copy1 - c^1_2_4_copy2 <= 0
 R101: d^1_4 - c^1_1_4_copy2 = 0
 R102: d^1_4 - c^1_2_4_copy2 = 0
 R103: e^1_4 - d^1_4 >= 0
 R104: e^1_4 - d^1_4 - c^1_1_3_copy1 = 0
 R105: - RK^0_5 + RK^0_5_copy1 <= 0
 R106: - RK^0_5 + RK^0_5_copy2 <= 0
 R107: RK^0_5 - RK^0_5_copy1 - RK^0_5_copy2 <= 0
 R108: - L^0_5 + L^0_5_copy1 <= 0
 R109: - L^0_5 + L^0_5_copy2 <= 0
 R110: L^0_5 - L^0_5_copy1 - L^0_5_copy2 <= 0
 R111: a^1_5 - RK^0_5_copy2 = 0
 R112: a^1_5 - L^0_5_copy2 = 0
 R113: - b^1_5 + R^0_5_copy2 = 0
 R114: - b^1_5 + c^1_1_4_copy1 >= 0
 R115: - a^1_5 - b^1_5 + c^1_1_4_copy1 = 0
 R116: - c^1_2_4_copy1 + c^1_1_5_copy1 <= 0
 R117: - c^1_2_4_copy1 + c^1_1_5_copy2 <= 0
 R118: c^1_2_4_copy1 - c^1_1_5_copy1 - c^1_1_5_copy2 <= 0
 R119: - c^1_7 + c^1_2_5_copy1 <= 0
 R120: - c^1_7 + c^1_2_5_copy2 <= 0
 R121: c^1_7 - c^1_2_5_copy1 - c^1_2_5_copy2 <= 0
 R122: d^1_5 - c^1_1_5_copy2 = 0
 R123: d^1_5 - c^1_2_5_copy2 = 0
 R124: e^1_5 - d^1_5 >= 0
 R125: e^1_5 - d^1_5 - c^1_1_4_copy1 = 0
 R126: - RK^0_6 + RK^0_6_copy1 <= 0
 R127: - RK^0_6 + RK^0_6_copy2 <= 0
 R128: RK^0_6 - RK^0_6_copy1 - RK^0_6_copy2 <= 0
 R129: - L^0_6 + L^0_6_copy1 <= 0
 R130: - L^0_6 + L^0_6_copy2 <= 0
 R131: L^0_6 - L^0_6_copy1 - L^0_6_copy2 <= 0
 R132: a^1_6 - RK^0_6_copy2 = 0
 R133: a^1_6 - L^0_6_copy2 = 0
 R134: - b^1_6 + R^0_6_copy2 = 0
 R135: - b^1_6 + c^1_1_5_copy1 >= 0
 R136: - a^1_6 - b^1_6 + c^1_1_5_copy1 = 0
 R137: - c^1_2_5_copy1 + c^1_1_6_copy1 <= 0
 R138: - c^1_2_5_copy1 + c^1_1_6_copy2 <= 0
 R139: c^1_2_5_copy1 - c^1_1_6_copy1 - c^1_1_6_copy2 <= 0
 R140: - c^1_8 + c^1_2_6_copy1 <= 0
 R141: - c^1_8 + c^1_2_6_copy2 <= 0
 R142: c^1_8 - c^1_2_6_copy1 - c^1_2_6_copy2 <= 0
 R143: d^1_6 - c^1_1_6_copy2 = 0
 R144: d^1_6 - c^1_2_6_copy2 = 0
 R145: e^1_6 - d^1_6 >= 0
 R146: e^1_6 - d^1_6 - c^1_1_5_copy1 = 0
 R147: - RK^0_7 + RK^0_7_copy1 <= 0
 R148: - RK^0_7 + RK^0_7_copy2 <= 0
 R149: RK^0_7 - RK^0_7_copy1 - RK^0_7_copy2 <= 0
 R150: - L^0_7 + L^0_7_copy1 <= 0
 R151: - L^0_7 + L^0_7_copy2 <= 0
 R152: L^0_7 - L^0_7_copy1 - L^0_7_copy2 <= 0
 R153: a^1_7 - RK^0_7_copy2 = 0
 R154: a^1_7 - L^0_7_copy2 = 0
 R155: - b^1_7 + R^0_7_copy2 = 0
 R156: - b^1_7 + c^1_1_6_copy1 >= 0
 R157: - a^1_7 - b^1_7 + c^1_1_6_copy1 = 0
 R158: - c^1_2_6_copy1 + c^1_1_7_copy1 <= 0
 R159: - c^1_2_6_copy1 + c^1_1_7_copy2 <= 0
 R160: c^1_2_6_copy1 - c^1_1_7_copy1 - c^1_1_7_copy2 <= 0
 R161: - c^1_9 + c^1_2_7_copy1 <= 0
 R162: - c^1_9 + c^1_2_7_copy2 <= 0
 R163: c^1_9 - c^1_2_7_copy1 - c^1_2_7_copy2 <= 0
 R164: d^1_7 - c^1_1_7_copy2 = 0
 R165: d^1_7 - c^1_2_7_copy2 = 0
 R166: e^1_7 - d^1_7 >= 0
 R167: e^1_7 - d^1_7 - c^1_1_6_copy1 = 0
 R168: - RK^0_8 + RK^0_8_copy1 <= 0
 R169: - RK^0_8 + RK^0_8_copy2 <= 0
 R170: RK^0_8 - RK^0_8_copy1 - RK^0_8_copy2 <= 0
 R171: - L^0_8 + L^0_8_copy1 <= 0
 R172: - L^0_8 + L^0_8_copy2 <= 0
 R173: L^0_8 - L^0_8_copy1 - L^0_8_copy2 <= 0
 R174: a^1_8 - RK^0_8_copy2 = 0
 R175: a^1_8 - L^0_8_copy2 = 0
 R176: - b^1_8 + R^0_8_copy2 = 0
 R177: - b^1_8 + c^1_1_7_copy1 >= 0
 R178: - a^1_8 - b^1_8 + c^1_1_7_copy1 = 0
 R179: - c^1_2_7_copy1 + c^1_1_8_copy1 <= 0
 R180: - c^1_2_7_copy1 + c^1_1_8_copy2 <= 0
 R181: c^1_2_7_copy1 - c^1_1_8_copy1 - c^1_1_8_copy2 <= 0
 R182: - c^1_10 + c^1_2_8_copy1 <= 0
 R183: - c^1_10 + c^1_2_8_copy2 <= 0
 R184: c^1_10 - c^1_2_8_copy1 - c^1_2_8_copy2 <= 0
 R185: d^1_8 - c^1_1_8_copy2 = 0
 R186: d^1_8 - c^1_2_8_copy2 = 0
 R187: e^1_8 - d^1_8 >= 0
 R188: e^1_8 - d^1_8 - c^1_1_7_copy1 = 0
 R189: - RK^0_9 + RK^0_9_copy1 <= 0
 R190: - RK^0_9 + RK^0_9_copy2 <= 0
 R191: RK^0_9 - RK^0_9_copy1 - RK^0_9_copy2 <= 0
 R192: - L^0_9 + L^0_9_copy1 <= 0
 R193: - L^0_9 + L^0_9_copy2 <= 0
 R194: L^0_9 - L^0_9_copy1 - L^0_9_copy2 <= 0
 R195: a^1_9 - RK^0_9_copy2 = 0
 R196: a^1_9 - L^0_9_copy2 = 0
 R197: - b^1_9 + R^0_9_copy2 = 0
 R198: - b^1_9 + c^1_1_8_copy1 >= 0
 R199: - a^1_9 - b^1_9 + c^1_1_8_copy1 = 0
 R200: - c^1_2_8_copy1 + c^1_1_9_copy1 <= 0
 R201: - c^1_2_8_copy1 + c^1_1_9_copy2 <= 0
 R202: c^1_2_8_copy1 - c^1_1_9_copy1 - c^1_1_9_copy2 <= 0
 R203: - c^1_11 + c^1_2_9_copy1 <= 0
 R204: - c^1_11 + c^1_2_9_copy2 <= 0
 R205: c^1_11 - c^1_2_9_copy1 - c^1_2_9_copy2 <= 0
 R206: d^1_9 - c^1_1_9_copy2 = 0
 R207: d^1_9 - c^1_2_9_copy2 = 0
 R208: e^1_9 - d^1_9 >= 0
 R209: e^1_9 - d^1_9 - c^1_1_8_copy1 = 0
 R210: - RK^0_10 + RK^0_10_copy1 <= 0
 R211: - RK^0_10 + RK^0_10_copy2 <= 0
 R212: RK^0_10 - RK^0_10_copy1 - RK^0_10_copy2 <= 0
 R213: - L^0_10 + L^0_10_copy1 <= 0
 R214: - L^0_10 + L^0_10_copy2 <= 0
 R215: L^0_10 - L^0_10_copy1 - L^0_10_copy2 <= 0
 R216: a^1_10 - RK^0_10_copy2 = 0
 R217: a^1_10 - L^0_10_copy2 = 0
 R218: - b^1_10 + R^0_10_copy2 = 0
 R219: - b^1_10 + c^1_1_9_copy1 >= 0
 R220: - a^1_10 - b^1_10 + c^1_1_9_copy1 = 0
 R221: - c^1_2_9_copy1 + c^1_1_10_copy1 <= 0
 R222: - c^1_2_9_copy1 + c^1_1_10_copy2 <= 0
 R223: c^1_2_9_copy1 - c^1_1_10_copy1 - c^1_1_10_copy2 <= 0
 R224: - c^1_12 + c^1_2_10_copy1 <= 0
 R225: - c^1_12 + c^1_2_10_copy2 <= 0
 R226: c^1_12 - c^1_2_10_copy1 - c^1_2_10_copy2 <= 0
 R227: d^1_10 - c^1_1_10_copy2 = 0
 R228: d^1_10 - c^1_2_10_copy2 = 0
 R229: e^1_10 - d^1_10 >= 0
 R230: e^1_10 - d^1_10 - c^1_1_9_copy1 = 0
 R231: - RK^0_11 + RK^0_11_copy1 <= 0
 R232: - RK^0_11 + RK^0_11_copy2 <= 0
 R233: RK^0_11 - RK^0_11_copy1 - RK^0_11_copy2 <= 0
 R234: - L^0_11 + L^0_11_copy1 <= 0
 R235: - L^0_11 + L^0_11_copy2 <= 0
 R236: L^0_11 - L^0_11_copy1 - L^0_11_copy2 <= 0
 R237: a^1_11 - RK^0_11_copy2 = 0
 R238: a^1_11 - L^0_11_copy2 = 0
 R239: - b^1_11 + R^0_11_copy2 = 0
 R240: - b^1_11 + c^1_1_10_copy1 >= 0
 R241: - a^1_11 - b^1_11 + c^1_1_10_copy1 = 0
 R242: - c^1_2_10_copy1 + c^1_1_11_copy1 <= 0
 R243: - c^1_2_10_copy1 + c^1_1_11_copy2 <= 0
 R244: c^1_2_10_copy1 - c^1_1_11_copy1 - c^1_1_11_copy2 <= 0
 R245: - c^1_13 + c^1_2_11_copy1 <= 0
 R246: - c^1_13 + c^1_2_11_copy2 <= 0
 R247: c^1_13 - c^1_2_11_copy1 - c^1_2_11_copy2 <= 0
 R248: d^1_11 - c^1_1_11_copy2 = 0
 R249: d^1_11 - c^1_2_11_copy2 = 0
 R250: e^1_11 - d^1_11 >= 0
 R251: e^1_11 - d^1_11 - c^1_1_10_copy1 = 0
 R252: - RK^0_12 + RK^0_12_copy1 <= 0
 R253: - RK^0_12 + RK^0_12_copy2 <= 0
 R254: RK^0_12 - RK^0_12_copy1 - RK^0_12_copy2 <= 0
 R255: - L^0_12 + L^0_12_copy1 <= 0
 R256: - L^0_12 + L^0_12_copy2 <= 0
 R257: L^0_12 - L^0_12_copy1 - L^0_12_copy2 <= 0
 R258: a^1_12 - RK^0_12_copy2 = 0
 R259: a^1_12 - L^0_12_copy2 = 0
 R260: - b^1_12 + R^0_12_copy2 = 0
 R261: - b^1_12 + c^1_1_11_copy1 >= 0
 R262: - a^1_12 - b^1_12 + c^1_1_11_copy1 = 0
 R263: - c^1_2_11_copy1 + c^1_1_12_copy1 <= 0
 R264: - c^1_2_11_copy1 + c^1_1_12_copy2 <= 0
 R265: c^1_2_11_copy1 - c^1_1_12_copy1 - c^1_1_12_copy2 <= 0
 R266: - c^1_14 + c^1_2_12_copy1 <= 0
 R267: - c^1_14 + c^1_2_12_copy2 <= 0
 R268: c^1_14 - c^1_2_12_copy1 - c^1_2_12_copy2 <= 0
 R269: d^1_12 - c^1_1_12_copy2 = 0
 R270: d^1_12 - c^1_2_12_copy2 = 0
 R271: e^1_12 - d^1_12 >= 0
 R272: e^1_12 - d^1_12 - c^1_1_11_copy1 = 0
 R273: - RK^0_13 + RK^0_13_copy1 <= 0
 R274: - RK^0_13 + RK^0_13_copy2 <= 0
 R275: RK^0_13 - RK^0_13_copy1 - RK^0_13_copy2 <= 0
 R276: - L^0_13 + L^0_13_copy1 <= 0
 R277: - L^0_13 + L^0_13_copy2 <= 0
 R278: L^0_13 - L^0_13_copy1 - L^0_13_copy2 <= 0
 R279: a^1_13 - RK^0_13_copy2 = 0
 R280: a^1_13 - L^0_13_copy2 = 0
 R281: - b^1_13 + R^0_13_copy2 = 0
 R282: - b^1_13 + c^1_1_12_copy1 >= 0
 R283: - a^1_13 - b^1_13 + c^1_1_12_copy1 = 0
 R284: - c^1_2_12_copy1 + c^1_1_13_copy1 <= 0
 R285: - c^1_2_12_copy1 + c^1_1_13_copy2 <= 0
 R286: c^1_2_12_copy1 - c^1_1_13_copy1 - c^1_1_13_copy2 <= 0
 R287: - c^1_15 + c^1_2_13_copy1 <= 0
 R288: - c^1_15 + c^1_2_13_copy2 <= 0
 R289: c^1_15 - c^1_2_13_copy1 - c^1_2_13_copy2 <= 0
 R290: d^1_13 - c^1_1_13_copy2 = 0
 R291: d^1_13 - c^1_2_13_copy2 = 0
 R292: e^1_13 - d^1_13 >= 0
 R293: e^1_13 - d^1_13 - c^1_1_12_copy1 = 0
 R294: - RK^0_14 + RK^0_14_copy1 <= 0
 R295: - RK^0_14 + RK^0_14_copy2 <= 0
 R296: RK^0_14 - RK^0_14_copy1 - RK^0_14_copy2 <= 0
 R297: - L^0_14 + L^0_14_copy1 <= 0
 R298: - L^0_14 + L^0_14_copy2 <= 0
 R299: L^0_14 - L^0_14_copy1 - L^0_14_copy2 <= 0
 R300: a^1_14 - RK^0_14_copy2 = 0
 R301: a^1_14 - L^0_14_copy2 = 0
 R302: - b^1_14 + R^0_14_copy2 = 0
 R303: - b^1_14 + c^1_1_13_copy1 >= 0
 R304: - a^1_14 - b^1_14 + c^1_1_13_copy1 = 0
 R305: e^1_14 - c^1_1_13_copy1 = 0
 R306: - RK^0_15 + RK^0_15_copy1 <= 0
 R307: - RK^0_15 + RK^0_15_copy2 <= 0
 R308: RK^0_15 - RK^0_15_copy1 - RK^0_15_copy2 <= 0
 R309: - L^0_15 + L^0_15_copy1 <= 0
 R310: - L^0_15 + L^0_15_copy2 <= 0
 R311: L^0_15 - L^0_15_copy1 - L^0_15_copy2 <= 0
 R312: a^1_15 - RK^0_15_copy2 = 0
 R313: a^1_15 - L^0_15_copy2 = 0
 R314: - b^1_15 + R^0_15_copy2 = 0
 R315: - b^1_15 + c^1_2_13_copy1 >= 0
 R316: - a^1_15 - b^1_15 + c^1_2_13_copy1 = 0
 R317: e^1_15 - c^1_2_13_copy1 = 0
 R318: - e^1_3 + e^1_3_0_copy1 <= 0
 R319: - e^1_3 + e^1_3_0_copy2 <= 0
 R320: e^1_3 - e^1_3_0_copy1 - e^1_3_0_copy2 <= 0
 R321: - e^1_9 + e^1_9_0_copy1 <= 0
 R322: - e^1_9 + e^1_9_0_copy2 <= 0
 R323: e^1_9 - e^1_9_0_copy1 - e^1_9_0_copy2 <= 0
 R324: - e^1_14 + e^1_14_0_copy1 <= 0
 R325: - e^1_14 + e^1_14_0_copy2 <= 0
 R326: e^1_14 - e^1_14_0_copy1 - e^1_14_0_copy2 <= 0
 R327: f^1_0 - e^1_3_0_copy2 - e^1_9_0_copy2 - e^1_14_0_copy2 = 0
 R328: f^1_0 - e^1_3_0_copy2 - e^1_9_0_copy2 >= 0
 R329: - e^1_4 + e^1_3_1_copy1 <= 0
 R330: - e^1_4 + e^1_3_1_copy2 <= 0
 R331: e^1_4 - e^1_3_1_copy1 - e^1_3_1_copy2 <= 0
 R332: - e^1_10 + e^1_9_1_copy1 <= 0
 R333: - e^1_10 + e^1_9_1_copy2 <= 0
 R334: e^1_10 - e^1_9_1_copy1 - e^1_9_1_copy2 <= 0
 R335: - e^1_15 + e^1_14_1_copy1 <= 0
 R336: - e^1_15 + e^1_14_1_copy2 <= 0
 R337: e^1_15 - e^1_14_1_copy1 - e^1_14_1_copy2 <= 0
 R338: f^1_1 - e^1_3_1_copy2 - e^1_9_1_copy2 - e^1_14_1_copy2 = 0
 R339: f^1_1 - e^1_3_1_copy2 - e^1_9_1_copy2 >= 0
 R340: - e^1_5 + e^1_3_2_copy1 <= 0
 R341: - e^1_5 + e^1_3_2_copy2 <= 0
 R342: e^1_5 - e^1_3_2_copy1 - e^1_3_2_copy2 <= 0
 R343: - e^1_11 + e^1_9_2_copy1 <= 0
 R344: - e^1_11 + e^1_9_2_copy2 <= 0
 R345: e^1_11 - e^1_9_2_copy1 - e^1_9_2_copy2 <= 0
 R346: - e^1_0 + e^1_14_2_copy1 <= 0
 R347: - e^1_0 + e^1_14_2_copy2 <= 0
 R348: e^1_0 - e^1_14_2_copy1 - e^1_14_2_copy2 <= 0
 R349: f^1_2 - e^1_3_2_copy2 - e^1_9_2_copy2 - e^1_14_2_copy2 = 0
 R350: f^1_2 - e^1_3_2_copy2 - e^1_9_2_copy2 >= 0
 R351: - e^1_6 + e^1_3_3_copy1 <= 0
 R352: - e^1_6 + e^1_3_3_copy2 <= 0
 R353: e^1_6 - e^1_3_3_copy1 - e^1_3_3_copy2 <= 0
 R354: - e^1_12 + e^1_9_3_copy1 <= 0
 R355: - e^1_12 + e^1_9_3_copy2 <= 0
 R356: e^1_12 - e^1_9_3_copy1 - e^1_9_3_copy2 <= 0
 R357: - e^1_1 + e^1_14_3_copy1 <= 0
 R358: - e^1_1 + e^1_14_3_copy2 <= 0
 R359: e^1_1 - e^1_14_3_copy1 - e^1_14_3_copy2 <= 0
 R360: f^1_3 - e^1_3_3_copy2 - e^1_9_3_copy2 - e^1_14_3_copy2 = 0
 R361: f^1_3 - e^1_3_3_copy2 - e^1_9_3_copy2 >= 0
 R362: - e^1_7 + e^1_3_4_copy1 <= 0
 R363: - e^1_7 + e^1_3_4_copy2 <= 0
 R364: e^1_7 - e^1_3_4_copy1 - e^1_3_4_copy2 <= 0
 R365: - e^1_13 + e^1_9_4_copy1 <= 0
 R366: - e^1_13 + e^1_9_4_copy2 <= 0
 R367: e^1_13 - e^1_9_4_copy1 - e^1_9_4_copy2 <= 0
 R368: - e^1_2 + e^1_14_4_copy1 <= 0
 R369: - e^1_2 + e^1_14_4_copy2 <= 0
 R370: e^1_2 - e^1_14_4_copy1 - e^1_14_4_copy2 <= 0
 R371: f^1_4 - e^1_3_4_copy2 - e^1_9_4_copy2 - e^1_14_4_copy2 = 0
 R372: f^1_4 - e^1_3_4_copy2 - e^1_9_4_copy2 >= 0
 R373: - e^1_8 + e^1_3_5_copy1 <= 0
 R374: - e^1_8 + e^1_3_5_copy2 <= 0
 R375: e^1_8 - e^1_3_5_copy1 - e^1_3_5_copy2 <= 0
 R376: - e^1_14_0_copy1 + e^1_9_5_copy1 <= 0
 R377: - e^1_14_0_copy1 + e^1_9_5_copy2 <= 0
 R378: e^1_14_0_copy1 - e^1_9_5_copy1 - e^1_9_5_copy2 <= 0
 R379: - e^1_3_0_copy1 + e^1_14_5_copy1 <= 0
 R380: - e^1_3_0_copy1 + e^1_14_5_copy2 <= 0
 R381: e^1_3_0_copy1 - e^1_14_5_copy1 - e^1_14_5_copy2 <= 0
 R382: f^1_5 - e^1_3_5_copy2 - e^1_9_5_copy2 - e^1_14_5_copy2 = 0
 R383: f^1_5 - e^1_3_5_copy2 - e^1_9_5_copy2 >= 0
 R384: - e^1_9_0_copy1 + e^1_3_6_copy1 <= 0
 R385: - e^1_9_0_copy1 + e^1_3_6_copy2 <= 0
 R386: e^1_9_0_copy1 - e^1_3_6_copy1 - e^1_3_6_copy2 <= 0
 R387: - e^1_14_1_copy1 + e^1_9_6_copy1 <= 0
 R388: - e^1_14_1_copy1 + e^1_9_6_copy2 <= 0
 R389: e^1_14_1_copy1 - e^1_9_6_copy1 - e^1_9_6_copy2 <= 0
 R390: - e^1_3_1_copy1 + e^1_14_6_copy1 <= 0
 R391: - e^1_3_1_copy1 + e^1_14_6_copy2 <= 0
 R392: e^1_3_1_copy1 - e^1_14_6_copy1 - e^1_14_6_copy2 <= 0
 R393: f^1_6 - e^1_3_6_copy2 - e^1_9_6_copy2 - e^1_14_6_copy2 = 0
 R394: f^1_6 - e^1_3_6_copy2 - e^1_9_6_copy2 >= 0
 R395: - e^1_9_1_copy1 + e^1_3_7_copy1 <= 0
 R396: - e^1_9_1_copy1 + e^1_3_7_copy2 <= 0
 R397: e^1_9_1_copy1 - e^1_3_7_copy1 - e^1_3_7_copy2 <= 0
 R398: - e^1_14_2_copy1 + e^1_9_7_copy1 <= 0
 R399: - e^1_14_2_copy1 + e^1_9_7_copy2 <= 0
 R400: e^1_14_2_copy1 - e^1_9_7_copy1 - e^1_9_7_copy2 <= 0
 R401: - e^1_3_2_copy1 + e^1_14_7_copy1 <= 0
 R402: - e^1_3_2_copy1 + e^1_14_7_copy2 <= 0
 R403: e^1_3_2_copy1 - e^1_14_7_copy1 - e^1_14_7_copy2 <= 0
 R404: f^1_7 - e^1_3_7_copy2 - e^1_9_7_copy2 - e^1_14_7_copy2 = 0
 R405: f^1_7 - e^1_3_7_copy2 - e^1_9_7_copy2 >= 0
 R406: - e^1_9_2_copy1 + e^1_3_8_copy1 <= 0
 R407: - e^1_9_2_copy1 + e^1_3_8_copy2 <= 0
 R408: e^1_9_2_copy1 - e^1_3_8_copy1 - e^1_3_8_copy2 <= 0
 R409: - e^1_14_3_copy1 + e^1_9_8_copy1 <= 0
 R410: - e^1_14_3_copy1 + e^1_9_8_copy2 <= 0
 R411: e^1_14_3_copy1 - e^1_9_8_copy1 - e^1_9_8_copy2 <= 0
 R412: - e^1_3_3_copy1 + e^1_14_8_copy1 <= 0
 R413: - e^1_3_3_copy1 + e^1_14_8_copy2 <= 0
 R414: e^1_3_3_copy1 - e^1_14_8_copy1 - e^1_14_8_copy2 <= 0
 R415: f^1_8 - e^1_3_8_copy2 - e^1_9_8_copy2 - e^1_14_8_copy2 = 0
 R416: f^1_8 - e^1_3_8_copy2 - e^1_9_8_copy2 >= 0
 R417: - e^1_9_3_copy1 + e^1_3_9_copy1 <= 0
 R418: - e^1_9_3_copy1 + e^1_3_9_copy2 <= 0
 R419: e^1_9_3_copy1 - e^1_3_9_copy1 - e^1_3_9_copy2 <= 0
 R420: - e^1_14_4_copy1 + e^1_9_9_copy1 <= 0
 R421: - e^1_14_4_copy1 + e^1_9_9_copy2 <= 0
 R422: e^1_14_4_copy1 - e^1_9_9_copy1 - e^1_9_9_copy2 <= 0
 R423: - e^1_3_4_copy1 + e^1_14_9_copy1 <= 0
 R424: - e^1_3_4_copy1 + e^1_14_9_copy2 <= 0
 R425: e^1_3_4_copy1 - e^1_14_9_copy1 - e^1_14_9_copy2 <= 0
 R426: f^1_9 - e^1_3_9_copy2 - e^1_9_9_copy2 - e^1_14_9_copy2 = 0
 R427: f^1_9 - e^1_3_9_copy2 - e^1_9_9_copy2 >= 0
 R428: - e^1_9_4_copy1 + e^1_3_10_copy1 <= 0
 R429: - e^1_9_4_copy1 + e^1_3_10_copy2 <= 0
 R430: e^1_9_4_copy1 - e^1_3_10_copy1 - e^1_3_10_copy2 <= 0
 R431: - e^1_14_5_copy1 + e^1_9_10_copy1 <= 0
 R432: - e^1_14_5_copy1 + e^1_9_10_copy2 <= 0
 R433: e^1_14_5_copy1 - e^1_9_10_copy1 - e^1_9_10_copy2 <= 0
 R434: - e^1_3_5_copy1 + e^1_14_10_copy1 <= 0
 R435: - e^1_3_5_copy1 + e^1_14_10_copy2 <= 0
 R436: e^1_3_5_copy1 - e^1_14_10_copy1 - e^1_14_10_copy2 <= 0
 R437: f^1_10 - e^1_3_10_copy2 - e^1_9_10_copy2 - e^1_14_10_copy2 = 0
 R438: f^1_10 - e^1_3_10_copy2 - e^1_9_10_copy2 >= 0
 R439: - e^1_9_5_copy1 + e^1_3_11_copy1 <= 0
 R440: - e^1_9_5_copy1 + e^1_3_11_copy2 <= 0
 R441: e^1_9_5_copy1 - e^1_3_11_copy1 - e^1_3_11_copy2 <= 0
 R442: - e^1_14_6_copy1 + e^1_9_11_copy1 <= 0
 R443: - e^1_14_6_copy1 + e^1_9_11_copy2 <= 0
 R444: e^1_14_6_copy1 - e^1_9_11_copy1 - e^1_9_11_copy2 <= 0
 R445: - e^1_3_6_copy1 + e^1_14_11_copy1 <= 0
 R446: - e^1_3_6_copy1 + e^1_14_11_copy2 <= 0
 R447: e^1_3_6_copy1 - e^1_14_11_copy1 - e^1_14_11_copy2 <= 0
 R448: f^1_11 - e^1_3_11_copy2 - e^1_9_11_copy2 - e^1_14_11_copy2 = 0
 R449: f^1_11 - e^1_3_11_copy2 - e^1_9_11_copy2 >= 0
 R450: - e^1_9_6_copy1 + e^1_3_12_copy1 <= 0
 R451: - e^1_9_6_copy1 + e^1_3_12_copy2 <= 0
 R452: e^1_9_6_copy1 - e^1_3_12_copy1 - e^1_3_12_copy2 <= 0
 R453: - e^1_14_7_copy1 + e^1_9_12_copy1 <= 0
 R454: - e^1_14_7_copy1 + e^1_9_12_copy2 <= 0
 R455: e^1_14_7_copy1 - e^1_9_12_copy1 - e^1_9_12_copy2 <= 0
 R456: - e^1_3_7_copy1 + e^1_14_12_copy1 <= 0
 R457: - e^1_3_7_copy1 + e^1_14_12_copy2 <= 0
 R458: e^1_3_7_copy1 - e^1_14_12_copy1 - e^1_14_12_copy2 <= 0
 R459: f^1_12 - e^1_3_12_copy2 - e^1_9_12_copy2 - e^1_14_12_copy2 = 0
 R460: f^1_12 - e^1_3_12_copy2 - e^1_9_12_copy2 >= 0
 R461: - e^1_9_7_copy1 + e^1_3_13_copy1 <= 0
 R462: - e^1_9_7_copy1 + e^1_3_13_copy2 <= 0
 R463: e^1_9_7_copy1 - e^1_3_13_copy1 - e^1_3_13_copy2 <= 0
 R464: - e^1_14_8_copy1 + e^1_9_13_copy1 <= 0
 R465: - e^1_14_8_copy1 + e^1_9_13_copy2 <= 0
 R466: e^1_14_8_copy1 - e^1_9_13_copy1 - e^1_9_13_copy2 <= 0
 R467: - e^1_3_8_copy1 + e^1_14_13_copy1 <= 0
 R468: - e^1_3_8_copy1 + e^1_14_13_copy2 <= 0
 R469: e^1_3_8_copy1 - e^1_14_13_copy1 - e^1_14_13_copy2 <= 0
 R470: f^1_13 - e^1_3_13_copy2 - e^1_9_13_copy2 - e^1_14_13_copy2 = 0
 R471: f^1_13 - e^1_3_13_copy2 - e^1_9_13_copy2 >= 0
 R472: - e^1_9_8_copy1 + e^1_3_14_copy1 <= 0
 R473: - e^1_9_8_copy1 + e^1_3_14_copy2 <= 0
 R474: e^1_9_8_copy1 - e^1_3_14_copy1 - e^1_3_14_copy2 <= 0
 R475: - e^1_14_9_copy1 + e^1_9_14_copy1 <= 0
 R476: - e^1_14_9_copy1 + e^1_9_14_copy2 <= 0
 R477: e^1_14_9_copy1 - e^1_9_14_copy1 - e^1_9_14_copy2 <= 0
 R478: - e^1_3_9_copy1 + e^1_14_14_copy1 <= 0
 R479: - e^1_3_9_copy1 + e^1_14_14_copy2 <= 0
 R480: e^1_3_9_copy1 - e^1_14_14_copy1 - e^1_14_14_copy2 <= 0
 R481: f^1_14 - e^1_3_14_copy2 - e^1_9_14_copy2 - e^1_14_14_copy2 = 0
 R482: f^1_14 - e^1_3_14_copy2 - e^1_9_14_copy2 >= 0
 R483: - e^1_9_9_copy1 + e^1_3_15_copy1 <= 0
 R484: - e^1_9_9_copy1 + e^1_3_15_copy2 <= 0
 R485: e^1_9_9_copy1 - e^1_3_15_copy1 - e^1_3_15_copy2 <= 0
 R486: - e^1_14_10_copy1 + e^1_9_15_copy1 <= 0
 R487: - e^1_14_10_copy1 + e^1_9_15_copy2 <= 0
 R488: e^1_14_10_copy1 - e^1_9_15_copy1 - e^1_9_15_copy2 <= 0
 R489: - e^1_3_10_copy1 + e^1_14_15_copy1 <= 0
 R490: - e^1_3_10_copy1 + e^1_14_15_copy2 <= 0
 R491: e^1_3_10_copy1 - e^1_14_15_copy1 - e^1_14_15_copy2 <= 0
 R492: f^1_15 - e^1_3_15_copy2 - e^1_9_15_copy2 - e^1_14_15_copy2 = 0
 R493: f^1_15 - e^1_3_15_copy2 - e^1_9_15_copy2 >= 0
 R494: - f^1_0 + f^1_0_copy1 <= 0
 R495: - f^1_0 + f^1_0_copy2 <= 0
 R496: f^1_0 - f^1_0_copy1 - f^1_0_copy2 <= 0
 R497: g^1_0 - f^1_0_copy2 = 0
 R498: g^1_0 - RK^0_0_copy1 = 0
 R499: L^1_0 - R^0_0_copy1 >= 0
 R500: - g^1_0 + L^1_0 - R^0_0_copy1 = 0
 R501: R^1_0 - L^0_0_copy1 >= 0
 R502: R^1_0 - L^0_0_copy1 - f^1_0_copy1 = 0
 R503: - f^1_1 + f^1_1_copy1 <= 0
 R504: - f^1_1 + f^1_1_copy2 <= 0
 R505: f^1_1 - f^1_1_copy1 - f^1_1_copy2 <= 0
 R506: g^1_1 - f^1_1_copy2 = 0
 R507: g^1_1 - RK^0_1_copy1 = 0
 R508: L^1_1 - R^0_1_copy1 >= 0
 R509: - g^1_1 + L^1_1 - R^0_1_copy1 = 0
 R510: R^1_1 - L^0_1_copy1 >= 0
 R511: R^1_1 - L^0_1_copy1 - f^1_1_copy1 = 0
 R512: - f^1_2 + f^1_2_copy1 <= 0
 R513: - f^1_2 + f^1_2_copy2 <= 0
 R514: f^1_2 - f^1_2_copy1 - f^1_2_copy2 <= 0
 R515: g^1_2 - f^1_2_copy2 = 0
 R516: g^1_2 - RK^0_2_copy1 = 0
 R517: L^1_2 - R^0_2_copy1 >= 0
 R518: - g^1_2 + L^1_2 - R^0_2_copy1 = 0
 R519: R^1_2 - L^0_2_copy1 >= 0
 R520: R^1_2 - L^0_2_copy1 - f^1_2_copy1 = 0
 R521: - f^1_3 + f^1_3_copy1 <= 0
 R522: - f^1_3 + f^1_3_copy2 <= 0
 R523: f^1_3 - f^1_3_copy1 - f^1_3_copy2 <= 0
 R524: g^1_3 - f^1_3_copy2 = 0
 R525: g^1_3 - RK^0_3_copy1 = 0
 R526: L^1_3 - R^0_3_copy1 >= 0
 R527: - g^1_3 + L^1_3 - R^0_3_copy1 = 0
 R528: R^1_3 - L^0_3_copy1 >= 0
 R529: R^1_3 - L^0_3_copy1 - f^1_3_copy1 = 0
 R530: - f^1_4 + f^1_4_copy1 <= 0
 R531: - f^1_4 + f^1_4_copy2 <= 0
 R532: f^1_4 - f^1_4_copy1 - f^1_4_copy2 <= 0
 R533: g^1_4 - f^1_4_copy2 = 0
 R534: g^1_4 - RK^0_4_copy1 = 0
 R535: L^1_4 - R^0_4_copy1 >= 0
 R536: - g^1_4 + L^1_4 - R^0_4_copy1 = 0
 R537: R^1_4 - L^0_4_copy1 >= 0
 R538: R^1_4 - L^0_4_copy1 - f^1_4_copy1 = 0
 R539: - f^1_5 + f^1_5_copy1 <= 0
 R540: - f^1_5 + f^1_5_copy2 <= 0
 R541: f^1_5 - f^1_5_copy1 - f^1_5_copy2 <= 0
 R542: g^1_5 - f^1_5_copy2 = 0
 R543: g^1_5 - RK^0_5_copy1 = 0
 R544: L^1_5 - R^0_5_copy1 >= 0
 R545: - g^1_5 + L^1_5 - R^0_5_copy1 = 0
 R546: R^1_5 - L^0_5_copy1 >= 0
 R547: R^1_5 - L^0_5_copy1 - f^1_5_copy1 = 0
 R548: - f^1_6 + f^1_6_copy1 <= 0
 R549: - f^1_6 + f^1_6_copy2 <= 0
 R550: f^1_6 - f^1_6_copy1 - f^1_6_copy2 <= 0
 R551: g^1_6 - f^1_6_copy2 = 0
 R552: g^1_6 - RK^0_6_copy1 = 0
 R553: L^1_6 - R^0_6_copy1 >= 0
 R554: - g^1_6 + L^1_6 - R^0_6_copy1 = 0
 R555: R^1_6 - L^0_6_copy1 >= 0
 R556: R^1_6 - L^0_6_copy1 - f^1_6_copy1 = 0
 R557: - f^1_7 + f^1_7_copy1 <= 0
 R558: - f^1_7 + f^1_7_copy2 <= 0
 R559: f^1_7 - f^1_7_copy1 - f^1_7_copy2 <= 0
 R560: g^1_7 - f^1_7_copy2 = 0
 R561: g^1_7 - RK^0_7_copy1 = 0
 R562: L^1_7 - R^0_7_copy1 >= 0
 R563: - g^1_7 + L^1_7 - R^0_7_copy1 = 0
 R564: R^1_7 - L^0_7_copy1 >= 0
 R565: R^1_7 - L^0_7_copy1 - f^1_7_copy1 = 0
 R566: - f^1_8 + f^1_8_copy1 <= 0
 R567: - f^1_8 + f^1_8_copy2 <= 0
 R568: f^1_8 - f^1_8_copy1 - f^1_8_copy2 <= 0
 R569: g^1_8 - f^1_8_copy2 = 0
 R570: g^1_8 - RK^0_8_copy1 = 0
 R571: L^1_8 - R^0_8_copy1 >= 0
 R572: - g^1_8 + L^1_8 - R^0_8_copy1 = 0
 R573: R^1_8 - L^0_8_copy1 >= 0
 R574: R^1_8 - L^0_8_copy1 - f^1_8_copy1 = 0
 R575: - f^1_9 + f^1_9_copy1 <= 0
 R576: - f^1_9 + f^1_9_copy2 <= 0
 R577: f^1_9 - f^1_9_copy1 - f^1_9_copy2 <= 0
 R578: g^1_9 - f^1_9_copy2 = 0
 R579: g^1_9 - RK^0_9_copy1 = 0
 R580: L^1_9 - R^0_9_copy1 >= 0
 R581: - g^1_9 + L^1_9 - R^0_9_copy1 = 0
 R582: R^1_9 - L^0_9_copy1 >= 0
 R583: R^1_9 - L^0_9_copy1 - f^1_9_copy1 = 0
 R584: - f^1_10 + f^1_10_copy1 <= 0
 R585: - f^1_10 + f^1_10_copy2 <= 0
 R586: f^1_10 - f^1_10_copy1 - f^1_10_copy2 <= 0
 R587: g^1_10 - f^1_10_copy2 = 0
 R588: g^1_10 - RK^0_10_copy1 = 0
 R589: L^1_10 - R^0_10_copy1 >= 0
 R590: - g^1_10 + L^1_10 - R^0_10_copy1 = 0
 R591: R^1_10 - L^0_10_copy1 >= 0
 R592: R^1_10 - L^0_10_copy1 - f^1_10_copy1 = 0
 R593: - f^1_11 + f^1_11_copy1 <= 0
 R594: - f^1_11 + f^1_11_copy2 <= 0
 R595: f^1_11 - f^1_11_copy1 - f^1_11_copy2 <= 0
 R596: g^1_11 - f^1_11_copy2 = 0
 R597: g^1_11 - RK^0_11_copy1 = 0
 R598: L^1_11 - R^0_11_copy1 >= 0
 R599: - g^1_11 + L^1_11 - R^0_11_copy1 = 0
 R600: R^1_11 - L^0_11_copy1 >= 0
 R601: R^1_11 - L^0_11_copy1 - f^1_11_copy1 = 0
 R602: - f^1_12 + f^1_12_copy1 <= 0
 R603: - f^1_12 + f^1_12_copy2 <= 0
 R604: f^1_12 - f^1_12_copy1 - f^1_12_copy2 <= 0
 R605: g^1_12 - f^1_12_copy2 = 0
 R606: g^1_12 - RK^0_12_copy1 = 0
 R607: L^1_12 - R^0_12_copy1 >= 0
 R608: - g^1_12 + L^1_12 - R^0_12_copy1 = 0
 R609: R^1_12 - L^0_12_copy1 >= 0
 R610: R^1_12 - L^0_12_copy1 - f^1_12_copy1 = 0
 R611: - f^1_13 + f^1_13_copy1 <= 0
 R612: - f^1_13 + f^1_13_copy2 <= 0
 R613: f^1_13 - f^1_13_copy1 - f^1_13_copy2 <= 0
 R614: g^1_13 - f^1_13_copy2 = 0
 R615: g^1_13 - RK^0_13_copy1 = 0
 R616: L^1_13 - R^0_13_copy1 >= 0
 R617: - g^1_13 + L^1_13 - R^0_13_copy1 = 0
 R618: R^1_13 - L^0_13_copy1 >= 0
 R619: R^1_13 - L^0_13_copy1 - f^1_13_copy1 = 0
 R620: - f^1_14 + f^1_14_copy1 <= 0
 R621: - f^1_14 + f^1_14_copy2 <= 0
 R622: f^1_14 - f^1_14_copy1 - f^1_14_copy2 <= 0
 R623: g^1_14 - f^1_14_copy2 = 0
 R624: g^1_14 - RK^0_14_copy1 = 0
 R625: L^1_14 - R^0_14_copy1 >= 0
 R626: - g^1_14 + L^1_14 - R^0_14_copy1 = 0
 R627: R^1_14 - L^0_14_copy1 >= 0
 R628: R^1_14 - L^0_14_copy1 - f^1_14_copy1 = 0
 R629: - f^1_15 + f^1_15_copy1 <= 0
 R630: - f^1_15 + f^1_15_copy2 <= 0
 R631: f^1_15 - f^1_15_copy1 - f^1_15_copy2 <= 0
 R632: g^1_15 - f^1_15_copy2 = 0
 R633: g^1_15 - RK^0_15_copy1 = 0
 R634: L^1_15 - R^0_15_copy1 >= 0
 R635: - g^1_15 + L^1_15 - R^0_15_copy1 = 0
 R636: R^1_15 - L^0_15_copy1 >= 0
 R637: R^1_15 - L^0_15_copy1 - f^1_15_copy1 = 0
 R638: - K_0 + RK^0_0 = 0
 R639: - flag_0 + RK^0_0 = 0
 R640: - K_1 + RK^0_1 = 0
 R641: - flag_1 + RK^0_1 = 0
 R642: - K_2 + RK^0_2 = 0
 R643: - flag_2 + RK^0_2 = 0
 R644: - K_3 + RK^0_3 = 0
 R645: - flag_3 + RK^0_3 = 0
 R646: - K_4 + RK^0_4 = 0
 R647: - flag_4 + RK^0_4 = 0
 R648: - K_5 + RK^0_5 = 0
 R649: - flag_5 + RK^0_5 = 0
 R650: - K_6 + RK^0_6 = 0
 R651: - flag_6 + RK^0_6 = 0
 R652: - K_7 + RK^0_7 = 0
 R653: - flag_7 + RK^0_7 = 0
 R654: - K_8 + RK^0_8 = 0
 R655: - flag_8 + RK^0_8 = 0
 R656: - K_9 + RK^0_9 = 0
 R657: - flag_9 + RK^0_9 = 0
 R658: - K_10 + RK^0_10 = 0
 R659: - flag_10 + RK^0_10 = 0
 R660: - K_11 + RK^0_11 = 0
 R661: - flag_11 + RK^0_11 = 0
 R662: - K_12 + RK^0_12 = 0
 R663: - flag_12 + RK^0_12 = 0
 R664: - K_13 + RK^0_13 = 0
 R665: - flag_13 + RK^0_13 = 0
 R666: - K_14 + RK^0_14 = 0
 R667: - flag_14 + RK^0_14 = 0
 R668: - K_15 + RK^0_15 = 0
 R669: - flag_15 + RK^0_15 = 0
 R670: L^0_0 = 1
 R671: L^0_1 = 1
 R672: L^0_2 = 1
 R673: L^0_3 = 1
 R674: L^0_4 = 1
 R675: L^0_5 = 1
 R676: L^0_6 = 1
 R677: L^0_7 = 1
 R678: L^0_8 = 1
 R679: L^0_9 = 1
 R680: L^0_10 = 1
 R681: L^0_11 = 1
 R682: L^0_12 = 1
 R683: L^0_13 = 1
 R684: L^0_14 = 1
 R685: L^0_15 = 1
 R686: R^0_0 = 1
 R687: R^0_1 = 1
 R688: R^0_2 = 1
 R689: R^0_3 = 1
 R690: R^0_4 = 1
 R691: R^0_5 = 1
 R692: R^0_6 = 1
 R693: R^0_7 = 1
 R694: R^0_8 = 1
 R695: R^0_9 = 1
 R696: R^0_10 = 1
 R697: R^0_11 = 1
 R698: R^0_12 = 1
 R699: R^0_13 = 1
 R700: R^0_14 = 1
 R701: R^0_15 = 0
 R702: L^1_0 = 1
 R703: R^1_0 = 0
 R704: L^1_1 = 0
 R705: R^1_1 = 0
 R706: L^1_2 = 0
 R707: R^1_2 = 0
 R708: L^1_3 = 0
 R709: R^1_3 = 0
 R710: L^1_4 = 0
 R711: R^1_4 = 0
 R712: L^1_5 = 0
 R713: R^1_5 = 0
 R714: L^1_6 = 0
 R715: R^1_6 = 0
 R716: L^1_7 = 0
 R717: R^1_7 = 0
 R718: L^1_8 = 0
 R719: R^1_8 = 0
 R720: L^1_9 = 0
 R721: R^1_9 = 0
 R722: L^1_10 = 0
 R723: R^1_10 = 0
 R724: L^1_11 = 0
 R725: R^1_11 = 0
 R726: L^1_12 = 0
 R727: R^1_12 = 0
 R728: L^1_13 = 0
 R729: R^1_13 = 0
 R730: L^1_14 = 0
 R731: R^1_14 = 0
 R732: L^1_15 = 0
 R733: R^1_15 = 0
Bounds
 flag_0 <= 2
 flag_1 <= 2
 flag_2 <= 2
 flag_3 <= 2
 flag_4 <= 2
 flag_5 <= 2
 flag_6 <= 2
 flag_7 <= 2
 flag_8 <= 2
 flag_9 <= 2
 flag_10 <= 2
 flag_11 <= 2
 flag_12 <= 2
 flag_13 <= 2
 flag_14 <= 2
 flag_15 <= 2
 flag_16 <= 2
 flag_17 <= 2
 flag_18 <= 2
 flag_19 <= 2
 flag_20 <= 2
 flag_21 <= 2
 flag_22 <= 2
 flag_23 <= 2
 flag_24 <= 2
 flag_25 <= 2
 flag_26 <= 2
 flag_27 <= 2
 flag_28 <= 2
 flag_29 <= 2
 flag_30 <= 2
 flag_31 <= 2
 flag_32 <= 2
 flag_33 <= 2
 flag_34 <= 2
 flag_35 <= 2
 flag_36 <= 2
 flag_37 <= 2
 flag_38 <= 2
 flag_39 <= 2
 flag_40 <= 2
 flag_41 <= 2
 flag_42 <= 2
 flag_43 <= 2
 flag_44 <= 2
 flag_45 <= 2
 flag_46 <= 2
 flag_47 <= 2
 flag_48 <= 2
 flag_49 <= 2
 flag_50 <= 2
 flag_51 <= 2
 flag_52 <= 2
 flag_53 <= 2
 flag_54 <= 2
 flag_55 <= 2
 flag_56 <= 2
 flag_57 <= 2
 flag_58 <= 2
 flag_59 <= 2
 flag_60 <= 2
 flag_61 <= 2
 flag_62 <= 2
 flag_63 <= 2
Binaries
 L^0_0 R^0_0 L^0_1 R^0_1 L^0_2 R^0_2 L^0_3 R^0_3 L^0_4 R^0_4 L^0_5 R^0_5
 L^0_6 R^0_6 L^0_7 R^0_7 L^0_8 R^0_8 L^0_9 R^0_9 L^0_10 R^0_10 L^0_11
 R^0_11 L^0_12 R^0_12 L^0_13 R^0_13 L^0_14 R^0_14 L^0_15 R^0_15 K_0 K_1 K_2
 K_3 K_4 K_5 K_6 K_7 K_8 K_9 K_10 K_11 K_12 K_13 K_14 K_15 K_16 K_17 K_18
 K_19 K_20 K_21 K_22 K_23 K_24 K_25 K_26 K_27 K_28 K_29 K_30 K_31 K_32 K_33
 K_34 K_35 K_36 K_37 K_38 K_39 K_40 K_41 K_42 K_43 K_44 K_45 K_46 K_47 K_48
 K_49 K_50 K_51 K_52 K_53 K_54 K_55 K_56 K_57 K_58 K_59 K_60 K_61 K_62 K_63
 RK^0_0 a^1_0 b^1_0 c^1_0 e^1_0 f^1_0 g^1_0 L^1_0 R^1_0 RK^0_1 a^1_1 b^1_1
 c^1_1 e^1_1 f^1_1 g^1_1 L^1_1 R^1_1 RK^0_2 a^1_2 b^1_2 c^1_2 e^1_2 f^1_2
 g^1_2 L^1_2 R^1_2 RK^0_3 a^1_3 b^1_3 c^1_3 e^1_3 f^1_3 g^1_3 L^1_3 R^1_3
 RK^0_4 a^1_4 b^1_4 c^1_4 e^1_4 f^1_4 g^1_4 L^1_4 R^1_4 RK^0_5 a^1_5 b^1_5
 c^1_5 e^1_5 f^1_5 g^1_5 L^1_5 R^1_5 RK^0_6 a^1_6 b^1_6 c^1_6 e^1_6 f^1_6
 g^1_6 L^1_6 R^1_6 RK^0_7 a^1_7 b^1_7 c^1_7 e^1_7 f^1_7 g^1_7 L^1_7 R^1_7
 RK^0_8 a^1_8 b^1_8 c^1_8 e^1_8 f^1_8 g^1_8 L^1_8 R^1_8 RK^0_9 a^1_9 b^1_9
 c^1_9 e^1_9 f^1_9 g^1_9 L^1_9 R^1_9 RK^0_10 a^1_10 b^1_10 c^1_10 e^1_10
 f^1_10 g^1_10 L^1_10 R^1_10 RK^0_11 a^1_11 b^1_11 c^1_11 e^1_11 f^1_11
 g^1_11 L^1_11 R^1_11 RK^0_12 a^1_12 b^1_12 c^1_12 e^1_12 f^1_12 g^1_12
 L^1_12 R^1_12 RK^0_13 a^1_13 b^1_13 c^1_13 e^1_13 f^1_13 g^1_13 L^1_13
 R^1_13 RK^0_14 a^1_14 b^1_14 c^1_14 e^1_14 f^1_14 g^1_14 L^1_14 R^1_14
 RK^0_15 a^1_15 b^1_15 c^1_15 e^1_15 f^1_15 g^1_15 L^1_15 R^1_15 d^1_0
 d^1_1 d^1_2 d^1_3 d^1_4 d^1_5 d^1_6 d^1_7 d^1_8 d^1_9 d^1_10 d^1_11 d^1_12
 d^1_13 RK^0_0_copy1 RK^0_0_copy2 L^0_0_copy1 L^0_0_copy2 R^0_0_copy1
 R^0_0_copy2 c^1_1_0_copy1 c^1_1_0_copy2 c^1_2_0_copy1 c^1_2_0_copy2
 RK^0_1_copy1 RK^0_1_copy2 L^0_1_copy1 L^0_1_copy2 R^0_1_copy1 R^0_1_copy2
 c^1_1_1_copy1 c^1_1_1_copy2 c^1_2_1_copy1 c^1_2_1_copy2 RK^0_2_copy1
 RK^0_2_copy2 L^0_2_copy1 L^0_2_copy2 R^0_2_copy1 R^0_2_copy2 c^1_1_2_copy1
 c^1_1_2_copy2 c^1_2_2_copy1 c^1_2_2_copy2 RK^0_3_copy1 RK^0_3_copy2
 L^0_3_copy1 L^0_3_copy2 R^0_3_copy1 R^0_3_copy2 c^1_1_3_copy1
 c^1_1_3_copy2 c^1_2_3_copy1 c^1_2_3_copy2 RK^0_4_copy1 RK^0_4_copy2
 L^0_4_copy1 L^0_4_copy2 R^0_4_copy1 R^0_4_copy2 c^1_1_4_copy1
 c^1_1_4_copy2 c^1_2_4_copy1 c^1_2_4_copy2 RK^0_5_copy1 RK^0_5_copy2
 L^0_5_copy1 L^0_5_copy2 R^0_5_copy1 R^0_5_copy2 c^1_1_5_copy1
 c^1_1_5_copy2 c^1_2_5_copy1 c^1_2_5_copy2 RK^0_6_copy1 RK^0_6_copy2
 L^0_6_copy1 L^0_6_copy2 R^0_6_copy1 R^0_6_copy2 c^1_1_6_copy1
 c^1_1_6_copy2 c^1_2_6_copy1 c^1_2_6_copy2 RK^0_7_copy1 RK^0_7_copy2
 L^0_7_copy1 L^0_7_copy2 R^0_7_copy1 R^0_7_copy2 c^1_1_7_copy1
 c^1_1_7_copy2 c^1_2_7_copy1 c^1_2_7_copy2 RK^0_8_copy1 RK^0_8_copy2
 L^0_8_copy1 L^0_8_copy2 R^0_8_copy1 R^0_8_copy2 c^1_1_8_copy1
 c^1_1_8_copy2 c^1_2_8_copy1 c^1_2_8_copy2 RK^0_9_copy1 RK^0_9_copy2
 L^0_9_copy1 L^0_9_copy2 R^0_9_copy1 R^0_9_copy2 c^1_1_9_copy1
 c^1_1_9_copy2 c^1_2_9_copy1 c^1_2_9_copy2 RK^0_10_copy1 RK^0_10_copy2
 L^0_10_copy1 L^0_10_copy2 R^0_10_copy1 R^0_10_copy2 c^1_1_10_copy1
 c^1_1_10_copy2 c^1_2_10_copy1 c^1_2_10_copy2 RK^0_11_copy1 RK^0_11_copy2
 L^0_11_copy1 L^0_11_copy2 R^0_11_copy1 R^0_11_copy2 c^1_1_11_copy1
 c^1_1_11_copy2 c^1_2_11_copy1 c^1_2_11_copy2 RK^0_12_copy1 RK^0_12_copy2
 L^0_12_copy1 L^0_12_copy2 R^0_12_copy1 R^0_12_copy2 c^1_1_12_copy1
 c^1_1_12_copy2 c^1_2_12_copy1 c^1_2_12_copy2 RK^0_13_copy1 RK^0_13_copy2
 L^0_13_copy1 L^0_13_copy2 R^0_13_copy1 R^0_13_copy2 c^1_1_13_copy1
 c^1_1_13_copy2 c^1_2_13_copy1 c^1_2_13_copy2 RK^0_14_copy1 RK^0_14_copy2
 L^0_14_copy1 L^0_14_copy2 R^0_14_copy1 R^0_14_copy2 RK^0_15_copy1
 RK^0_15_copy2 L^0_15_copy1 L^0_15_copy2 R^0_15_copy1 R^0_15_copy2
 e^1_3_0_copy1 e^1_3_0_copy2 e^1_9_0_copy1 e^1_9_0_copy2 e^1_14_0_copy1
 e^1_14_0_copy2 e^1_3_1_copy1 e^1_3_1_copy2 e^1_9_1_copy1 e^1_9_1_copy2
 e^1_14_1_copy1 e^1_14_1_copy2 e^1_3_2_copy1 e^1_3_2_copy2 e^1_9_2_copy1
 e^1_9_2_copy2 e^1_14_2_copy1 e^1_14_2_copy2 e^1_3_3_copy1 e^1_3_3_copy2
 e^1_9_3_copy1 e^1_9_3_copy2 e^1_14_3_copy1 e^1_14_3_copy2 e^1_3_4_copy1
 e^1_3_4_copy2 e^1_9_4_copy1 e^1_9_4_copy2 e^1_14_4_copy1 e^1_14_4_copy2
 e^1_3_5_copy1 e^1_3_5_copy2 e^1_9_5_copy1 e^1_9_5_copy2 e^1_14_5_copy1
 e^1_14_5_copy2 e^1_3_6_copy1 e^1_3_6_copy2 e^1_9_6_copy1 e^1_9_6_copy2
 e^1_14_6_copy1 e^1_14_6_copy2 e^1_3_7_copy1 e^1_3_7_copy2 e^1_9_7_copy1
 e^1_9_7_copy2 e^1_14_7_copy1 e^1_14_7_copy2 e^1_3_8_copy1 e^1_3_8_copy2
 e^1_9_8_copy1 e^1_9_8_copy2 e^1_14_8_copy1 e^1_14_8_copy2 e^1_3_9_copy1
 e^1_3_9_copy2 e^1_9_9_copy1 e^1_9_9_copy2 e^1_14_9_copy1 e^1_14_9_copy2
 e^1_3_10_copy1 e^1_3_10_copy2 e^1_9_10_copy1 e^1_9_10_copy2
 e^1_14_10_copy1 e^1_14_10_copy2 e^1_3_11_copy1 e^1_3_11_copy2
 e^1_9_11_copy1 e^1_9_11_copy2 e^1_14_11_copy1 e^1_14_11_copy2
 e^1_3_12_copy1 e^1_3_12_copy2 e^1_9_12_copy1 e^1_9_12_copy2
 e^1_14_12_copy1 e^1_14_12_copy2 e^1_3_13_copy1 e^1_3_13_copy2
 e^1_9_13_copy1 e^1_9_13_copy2 e^1_14_13_copy1 e^1_14_13_copy2
 e^1_3_14_copy1 e^1_3_14_copy2 e^1_9_14_copy1 e^1_9_14_copy2
 e^1_14_14_copy1 e^1_14_14_copy2 e^1_3_15_copy1 e^1_3_15_copy2
 e^1_9_15_copy1 e^1_9_15_copy2 e^1_14_15_copy1 e^1_14_15_copy2 f^1_0_copy1
 f^1_0_copy2 f^1_1_copy1 f^1_1_copy2 f^1_2_copy1 f^1_2_copy2 f^1_3_copy1
 f^1_3_copy2 f^1_4_copy1 f^1_4_copy2 f^1_5_copy1 f^1_5_copy2 f^1_6_copy1
 f^1_6_copy2 f^1_7_copy1 f^1_7_copy2 f^1_8_copy1 f^1_8_copy2 f^1_9_copy1
 f^1_9_copy2 f^1_10_copy1 f^1_10_copy2 f^1_11_copy1 f^1_11_copy2
 f^1_12_copy1 f^1_12_copy2 f^1_13_copy1 f^1_13_copy2 f^1_14_copy1
 f^1_14_copy2 f^1_15_copy1 f^1_15_copy2
Generals
 flag_0 flag_1 flag_2 flag_3 flag_4 flag_5 flag_6 flag_7 flag_8 flag_9
 flag_10 flag_11 flag_12 flag_13 flag_14 flag_15 flag_16 flag_17 flag_18
 flag_19 flag_20 flag_21 flag_22 flag_23 flag_24 flag_25 flag_26 flag_27
 flag_28 flag_29 flag_30 flag_31 flag_32 flag_33 flag_34 flag_35 flag_36
 flag_37 flag_38 flag_39 flag_40 flag_41 flag_42 flag_43 flag_44 flag_45
 flag_46 flag_47 flag_48 flag_49 flag_50 flag_51 flag_52 flag_53 flag_54
 flag_55 flag_56 flag_57 flag_58 flag_59 flag_60 flag_61 flag_62 flag_63
End
