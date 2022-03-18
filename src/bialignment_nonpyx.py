from collections import defaultdict

__version__ = "0.3b5"

blosum62 = """-  A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4
B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4
Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4
* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1
"""


def read_simmatrix(filename, scale=100):
    if filename == "BLOSUM62":
        lines = blosum62.split("\n")
    else:
        with open(filename, "r") as fh:
            lines = fh.readlines()

    keys = None
    keys2 = []
    matrix = dict()

    for i, line in enumerate(lines):
        if keys and i > len(keys):
            break
        line = line.split()
        if line[0] == "-":
            keys = line[1:]
        else:
            keys2.append(line[0])
            matrix[line[0]] = {
                key: (scale * int(val)) for key, val in zip(keys, line[1:])
            }

    if not keys == keys2:
        print("ERROR while reading simmatrix {filename}.")
    return matrix


def read_molecule(content, type):
    if type != "Protein":
        raise IOError(f"Cannot read files of type {type}")

    result = defaultdict(lambda: "")
    keys = ["Query", "Struc"]
    for line in content.split("\n"):
        line = line.split()
        if not line:
            continue
        if line[0] in keys:
            if len(line) != 4:
                raise IOError("Cannot parse")
            result[line[0]] += line[2]

    if len(result[keys[0]]) != len(result[keys[1]]):
        raise IOError("Sequence and structure of unequal length.")
    if len(result[keys[0]]) == 0:
        raise IOError("Input does not contain input sequence and structure.")

    return [result[k] for k in keys]


def read_molecule_from_file(filename, type):
    try:
        with open(filename, "r") as fh:
            return read_molecule(fh.read(), type)
    except FileNotFoundError as e:
        print("Input file not found.")
        print(e)
        sys.exit(-1)
    except IOError as e:
        print(f"Cannot read input file {filename}.")
        print(e)
        sys.exit(-1)


def breaklines(alilines, width):
    res = []

    offset = 0
    length = len(alilines[0][1])

    while offset < length:
        resblock = []
        for i, line in enumerate(alilines):
            name = line[0]
            aliline = line[1][offset : offset + width]
            resblock.append((name, aliline))
        offset += width

        res.append(resblock)

    return res


def runs(s):
    if s == "":
        return

    last_start = 0
    last = s[0]
    for i, x in enumerate(s[1:]):
        if x != last:
            yield (last, last_start, i + 1)
            last_start = i + 1
            last = x
    yield (last, last_start, len(s))


# struct = "HHHCCC---TTTCCHHHHHTT--EEECC"
# print(len(struct))
# for c,s,e in runs(struct):
#    print(c,s,e,struct[s:e])

helix_yadd_a = []
helix_yadd_b = []


def fourway_from_full(alilines):
    return [alilines[i] for i in [1, 3, 6, 8, 12, 13]]


def plot_alignment(
    alilines,
    width,
    *,
    show_structure_strings=False,
    name_offset=12,
    show_position_numbers=True,
    show_inconcruence = True,
    outname=None,
):
    """
    Plot bialignment

    Plots the bialignment with in several ways, optionally write plot to file

    Args:
        alilines : alignment strings with names; [(nameA, seqA), (nameB, seqB), ('',strA), ('',strB)]
        width : width of alignment, if longer break into rows
        show_structure_strings : switch display of structure strings
        name_offset : offset for showing names; should be adapted to name length
        show_positions_numbers : switch display of numbers at start and end of rows
        show_incongruence : switch display of incongruency bars
        outname : name of output file (extension controls format; see matplotlib.pyplot.savefig)
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    from collections import defaultdict

    helix_yadd_a = [0.0075]
    helix_yadd_b = [0.0075]

    if len(alilines) >= 13:  # full alilines
        alilines = fourway_from_full(alilines)

    aliblocks = breaklines(alilines, width)

    numblocks = len(aliblocks)
    fig, axs = plt.subplots(numblocks, 1, figsize=(0.18 * width, 2 * numblocks))

    if numblocks == 1:
        axs = [axs]

    font = {"family": "monospace", "weight": "normal", "size": 16.0}
    plt.rc("font", **font)  # pass in the font dict as kwargs

    colors = defaultdict(lambda: "grey", E="green", C="orange", T="blue", H="red")
    colors["-"] = None

    def draw_line(ax, s, e, y, color, lw):
        ax.plot(
            [s, e],
            [y + 0.025, y + 0.025],
            linewidth=lw,
            color=color,
            solid_capstyle="butt",
        )
        # lwcorr = (8-lw)/20
        # ax.plot([s+0.48-lwcorr,e-0.36+lwcorr], [y+0.025,y+0.025], linewidth=lw, color=color)

    def draw_sheet(ax, s, e, y, color):
        if s + 1 < e:
            ax.plot(
                [s, e - 1],
                [y + 0.025, y + 0.025],
                linewidth=8,
                color=color,
                solid_capstyle="butt",
            )
        ax.plot(
            [e - 0.05], [y + 0.025], linewidth=0, color=color, marker=5, markersize=13
        )

    def draw_helix(ax, s, e, y, color, yadd):
        xs = list(reversed(range(s, e + 1)))
        y += 0.025
        ys = [y + yadd[0]]
        for i in reversed(range(s, e)):
            yadd[0] = -yadd[0]
            ys.append(y + yadd[0])

        ax.plot(
            xs,
            ys,
            linewidth=6,
            color=color,
            solid_capstyle="butt",
            solid_joinstyle="round",
        )

    def draw_str(ax, y, line, yadd):
        name, struc = line
        for ch, s, e in reversed(list(runs(struc))):
            color = colors[ch]
            if ch == "E":
                draw_sheet(ax, s, e, y, color)
            elif ch == "H":
                draw_helix(ax, s, e, y, color, yadd)
            elif color == None:
                pass
            else:
                if ch == "T":
                    w = 8
                else:
                    w = 4
                draw_line(ax, s, e, y, color, w)

    def draw_seq(ax, y, line, other=None):
        name, seq = line
        ax.text(-name_offset, y, name)
        for x, ch in enumerate(seq):
            weight = "normal"
            color = "black"
            if other is not None:
                if ch != "-" and other[x] != "-":
                    color = "darkred"
                if ch == other[x]:
                    weight = "bold"
                    color = "black"
            ax.text(x, y, ch, weight=weight, color=color)

    def draw_shifts(ax, aa, bb):
        for x, (a, b) in enumerate(zip(aa, bb)):
            if a in ["<", ">"] or b in ["<", ">"]:
                ax.add_patch(
                    Rectangle(
                        (x, -0.022), 1, 0.4, edgecolor="black", fill=False, lw=0.5
                    )
                )

    ## store info on incongruence in former blocks/alignment rows
    incongruence_info = [0,0]

    def draw_incongruence(ax, aa, bb):

        def draw_incongruence_single(k, s, e, num):
            y = -0.0425 if k==1 else 0.405

            num += 0 # debug
            if num == 0:
                return

            color = 'darkred'
            if num<0:
                color = 'darkblue'

            num = abs(num)

            if s > e:
                return

            for i in range(num):
                o=0
                if num>1:
                    o = (i/(num-1)-0.5) * 0.02
                ax.plot(
                    [s, e+1],
                    [y+o, y+o],
                    linewidth=1,
                    color=color,
                    solid_capstyle="butt",
                )

        starts = [0,0]

        for x, ab in enumerate(zip(aa, bb)):
            for k,c in enumerate(ab):
                if c in ['<','>']:
                    draw_incongruence_single(k,starts[k],x-1,incongruence_info[k])
                    starts[k] = x+1
                    incongruence_info[k] += 1 if c=='>' else -1

        for k in range(2):
            draw_incongruence_single(k,starts[k],x,incongruence_info[k])


    offset_a = 1
    offset_b = 1

    for k, block in enumerate(aliblocks):
        ax = axs[k]
        ax.set_xlim(-0.5, width + 0.5)
        ax.set_ylim(-0.175, 0.425)
        ax.axis("off")

        # print(k, block)

        # for i, (name, aliline) in enumerate(block):
        #    print(f"{i:2} {name:12} {aliline}")
        # print()

        length = len(block[0][1])
        length_a = len(block[0][1].replace("-", ""))
        length_b = len(block[1][1].replace("-", ""))

        ## print sequence positions
        if show_position_numbers:
            ax.text(0, 0.435, offset_a, fontsize=10)
            offset_a += length_a
            ax.text(length, 0.435, offset_a - 1, fontsize=10, ha="right")

            ax.text(0, -0.12, offset_b, fontsize=10)
            offset_b += length_b
            ax.text(length, -0.12, offset_b - 1, fontsize=10, ha="right")

        draw_seq(ax, 0.2, block[0], block[1][1])
        draw_seq(ax, 0.1, block[1], block[0][1])

        draw_str(ax, 0.3, block[2], helix_yadd_a)
        draw_str(ax, 0.025, block[3], helix_yadd_b)
        if show_structure_strings:
            draw_seq(ax, 0.3, ("", block[2][1]))
            draw_seq(ax, 0, ("", block[3][1]))

    if len(block) > 4:
        for k, block in enumerate(aliblocks):
            ax = axs[k]
            draw_seq(ax, 0.375, ("", block[4][1].replace(".", " ")))
            draw_seq(ax, -0.075, ("", block[5][1].replace(".", " ")))
            draw_shifts(ax, block[4][1], block[5][1])
            draw_incongruence(ax, block[4][1], block[5][1])

    if outname is not None:
        plt.savefig(outname)
    plt.show()
