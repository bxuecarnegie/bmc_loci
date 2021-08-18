from graphics import *


def draw_operon(type_assignment, gene_list, organism, locusname, total_protein_length):
    # draw operon on screen using graphics library
    # protein type vs color definitions
    # First item: if pfam field empty
    # Second item: alc_dh
    # Third item: alddh
    # Forth item: signature enzymes
    # Sixth item: ptac

    type_list = [[("-",), "grey"], [("alcdh",), "green"], [("HMMalddh",), "red"],
                 [("HMMpropdeh", "HMMetly", "HMMgre", "HMMaldol", "HMMspu"), "purple"], [("HMMreg",), "orange"],
                 [("HMMptac", "HMM01515PTA_B"), "magenta"], [("H_", "Hp_"), "blue"],
                 [("T_", "Tsp_", "Ts_"), "cyan"], [("Tdp_",), (color_rgb(0, 160, 160))], [("P_",), "yellow"]]

    # read in operon
    if total_protein_length < 3000:
        total_protein_length = 3000  # for scaling
    width = 1400
    height = 700
    win = GraphWin(organism, width + 250, height, autoflush=False)  # give title and dimensions
    win.setBackground('white')
    win.update()
    current_pos = 0
    offset = 52
    previous_gene_number = int(gene_list[0][0])
    for field in gene_list:
        id_number = field[0]
        gene_length = int(width * (field[1] / total_protein_length))  # scaled to overall width
        assignment = field[2]
        hmm_bmc_type = field[3]
        direction = field[4]
        ep_string = field[5]
        id_string = id_number + " " + assignment + " " + hmm_bmc_type + " " + ep_string
        id_string_clean = id_string.replace("()", "").replace("(-)", "")
        difference = int(id_number) - previous_gene_number
        if abs(difference) > 50:  # if large break in gene id numbers, insert spacer
            offset = offset + 10
            diff_label = Text(Point(current_pos, offset), "(.." + str(difference) + "..)")
            diff_label.draw(win)
            offset = offset + 26
        # drawing of bar/arrow for single gene
        point1 = Point(current_pos, offset)
        point2 = Point(current_pos + gene_length, offset)
        rect = Line(point1, point2)
        rect.setWidth(15)
        if direction == "-":  # directionality of arrows
            rect.setArrow('first')
        if direction == "+":
            rect.setArrow('last')
        if (direction != "-") and (direction != "+"):
            rect.setArrow('none')
        bar_color = 'black'
        for pftype in type_list:
            for pftypemember in pftype[0]:
                if pftypemember in assignment:
                    bar_color = pftype[1]
        rect.setFill(bar_color)
        rect.draw(win)
        # drawing circle for encapsulation peptide position
        center_pos = Point(current_pos + int(gene_length / 2), offset)
        left_pos = Point(current_pos, offset)
        right_pos = Point(current_pos + int(gene_length), offset)
        if ep_string != "":
            ep_pos = center_pos
            if "NTERM" in ep_string:
                if direction == "-":
                    ep_pos = right_pos
                else:
                    ep_pos = left_pos
            if "CTERM" in ep_string:
                if direction == "-":
                    ep_pos = left_pos
                else:
                    ep_pos = right_pos
            ep_circle = Circle(ep_pos, 6)
            ep_circle.setFill('brown')
            ep_circle.setOutline('yellow')
            ep_circle.draw(win)
        label = Text(Point(current_pos + gene_length + len(id_string_clean) * 5, offset), id_string_clean)
        # coloring text green if type assignment matches best type from HMM match
        if type_assignment in hmm_bmc_type:
            label.setFill('green')
        else:
            if "-" in field[3]:
                label.setFill('grey')
            else:
                label.setFill('red')
        label.draw(win)
        # calculate difference between gene numbers and show that on right side
        previous_gene_number = int(id_number)
        current_pos = current_pos + gene_length
        offset = offset + 16
        if offset > height - 60:  # new column if at bottom of screen
            offset = 100
    # drawing legend
    legend_start = 250
    ls = legend_start
    message1 = Text(Point(89, ls), "BMC-P")
    message1.draw(win)
    rect = Line(Point(40, ls - 10), Point(40, ls + 10))
    rect.setWidth(18)
    rect.setFill('yellow')
    rect.draw(win)

    ls = legend_start + 20
    message1 = Text(Point(89, ls), "BMC-H")
    message1.draw(win)
    rect = Line(Point(40, ls - 10), Point(40, ls + 10))
    rect.setWidth(18)
    rect.setFill('blue')
    rect.draw(win)

    ls = legend_start + 40
    message1 = Text(Point(89, ls), "BMC-T")
    message1.draw(win)
    rect = Line(Point(40, ls - 10), Point(40, ls + 10))
    rect.setWidth(18)
    rect.setFill('cyan')
    rect.draw(win)

    ls = legend_start + 60
    message1 = Text(Point(155, ls), "Aldehyde dehydrogenase")
    message1.draw(win)
    rect = Line(Point(40, ls - 10), Point(40, ls + 10))
    rect.setWidth(18)
    rect.setFill('red')
    rect.draw(win)

    ls = legend_start + 80
    message1 = Text(Point(87, ls), "PTAC")
    message1.draw(win)
    rect = Line(Point(40, ls - 10), Point(40, ls + 10))
    rect.setWidth(18)
    rect.setFill('magenta')
    rect.draw(win)

    ls = legend_start + 100
    message1 = Text(Point(129, ls), "Signature enzyme")
    message1.draw(win)
    rect = Line(Point(40, ls - 10), Point(40, ls + 10))
    rect.setWidth(18)
    rect.setFill('purple')
    rect.draw(win)

    ls = legend_start + 120
    message1 = Text(Point(149, ls), "Alcohol dehydrogenase")
    message1.draw(win)
    rect = Line(Point(40, ls - 10), Point(40, ls + 10))
    rect.setWidth(18)
    rect.setFill('green')
    rect.draw(win)

    ls = legend_start + 140
    message1 = Text(Point(100, ls), "Regulator")
    message1.draw(win)
    rect = Line(Point(40, ls - 10), Point(40, ls + 10))
    rect.setWidth(18)
    rect.setFill('orange')
    rect.draw(win)

    message1 = Text(Point((width + 250) / 2, 20), organism)
    message2 = Text(Point((width + 250) / 2, 60), type_assignment)
    message3 = Text(Point((width + 250) / 2, 40), locusname)
    message1.draw(win)
    message2.setFill('green')
    message2.draw(win)
    message3.draw(win)
    # win.postscript(file=locusname+".eps")  output to eps file
    win.getKey()
    win.close()
