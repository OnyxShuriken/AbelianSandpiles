from tkinter import *
from tkinter import ttk
import math

for abc in range(1,28):
    # Create the main window
    root = Tk()
    root.title("Binary Hyperbolic " + str(abc))
    # Set window dimensions
    root.geometry("800x800")

    # Set the initial canvas size
    canvas_width = 805
    canvas_height = 805

    canvas = Canvas(root, width=canvas_width, height=canvas_height, 
    bg="white")
    canvas.pack()

    data = open("binaryHyperbolic/" + str(abc) + ".txt").read().split(',')

    colours = {'0':"grey1", '1':"grey25", '2':"grey50", '3':"grey75", 
    '4':"grey100"}

    length = int(math.log(len(data) + 1, 2))

    box_height = canvas_height / length

    y_offset = 0

    for i in range(0, length):
        x_offset = 0
        for ii in range(2**i, 2**(i+1)):
            real = ii - 1        
            canvas.create_rectangle(x_offset, y_offset,
                                    x_offset + (canvas_width / (2**(i))), 
    y_offset + box_height,
                                    fill=colours[data[real]], 
    outline=colours[data[real]])
            x_offset += canvas_width / (2**(i))
        y_offset += box_height


    # Run the Tkinter event loop
    root.mainloop() 
