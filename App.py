from tkinter import *
import tkinter as tk
from Algorithms import brute_force_search
from Algorithms import knuth_morris_pratt
from Algorithms import rabin_karp
import os.path
from tkinter import messagebox


"""This is the path variable.
   The user will need to change this variable depending on where your file
   and it also depends on your operating system
   
   Important:  MAKE SURE THAT THE FILE YOU WANT TO TEST ON IS PRESENT IN THE Genomes folder
                
               If it is not present, then set PATH = "" , and Specify the complete directory where your file exists 
"""


PATH = "Genomes\\"


class ProjectApp:
    def __init__(self):
        self.screen = None
        self.algorithms = {0: brute_force_search, 1: knuth_morris_pratt, 2: rabin_karp}
        self.width = 1000
        self.height = 400
        self.root = Tk()
        self.root.geometry(str(self.width) + "x" + str(self.height))
        self.light_text = "Enter your Query"

        self.root.title("Algorithms (CZ2001) Project 1")

        """These are the files which we gave the user as options.
            The user can also upload a file of their choice if they wish to test it on other files
        """

        self.files = {0: "GCA_000853505.1_ViralProj14913_genomic.fna",
                      1: "GCA_900205325.1_CV777_WBR_genomic.fna",
                      2: "GCA_011537355.1_ASM1153735v1_genomic.fna",
                      3: "GCA_900197115.1_PEDV_GER_L00857-K14_14-04_2014_genomic.fna",
                      4: "GCA_001904885.1_ViralProj357487_genomic.fna"}

        self.filename = ""
        self.genome_sequence = ""
        self.query = ""
        self.sequence_length = 0
        self.query_length = 0

        # variables to store user's choices
        self.file_choice = tk.IntVar()
        self.algorithm_choice = tk.IntVar()

        self.menu_bar = tk.Menu(self.root)

        self.menu_algorithm = tk.Menu(self.root)
        self.menu_algorithm.add_radiobutton(label="Brute Force", command=lambda: self.update_algorithm_choice(0))
        self.menu_algorithm.add_radiobutton(label="KMP", command=lambda: self.update_algorithm_choice(1))
        self.menu_algorithm.add_radiobutton(label="Rabin Karp", command=lambda: self.update_algorithm_choice(2))

        self.menu_file = tk.Menu(self.root)
        self.menu_file.add_radiobutton(label="File 1", command=lambda: self.update_file_choice(0))
        self.menu_file.add_radiobutton(label="File 2", command=lambda: self.update_file_choice(1))
        self.menu_file.add_radiobutton(label="File 3", command=lambda: self.update_file_choice(2))
        self.menu_file.add_radiobutton(label="File 4", command=lambda: self.update_file_choice(3))
        self.menu_file.add_radiobutton(label="File 5", command=lambda: self.update_file_choice(4))
        self.menu_file.add_radiobutton(label="Upload File", command=self.upload_file)

        self.menu_bar.add_cascade(label="File", menu=self.menu_file)
        self.menu_bar.add_cascade(label="Algorithm", menu=self.menu_algorithm)

        self.clear_button = Button(self.root, text="Clear", height=1, width=8, bg="RED", command=self.clear_text,
                                   activebackground='#e34234')
        self.search_button = Button(self.root, text="Search", height=1, width=8, bg="GREEN", activebackground='#44ff44',
                                    command=self.perform_search)

        self.scroll_bar = Scrollbar(self.root)
        self.entry_field = Entry(self.root, width=60, font=("Times New Roman", 20, "bold"), justify="center",
                                 bg="WHITE")
        self.entry_field.insert(0, self.light_text)
        self.entry_field.config(fg='grey')
        self.entry_field.bind('<FocusIn>', self.on_entry_click)
        self.entry_field.bind('<FocusOut>', self.on_focusout)
        self.result_output = Text(self.root, height=16, width=100)

        self.place_widgets()

        self.root.resizable(0, 0)
        self.root.mainloop()

    def place_widgets(self):
        self.root.config(menu=self.menu_bar)
        self.entry_field.place(x=100, y=10)
        self.clear_button.place(x=876, y=60)
        self.search_button.place(x=100, y=60)
        self.scroll_bar.config(command=self.result_output.yview)
        self.result_output.config(yscrollcommand=self.scroll_bar.set)

        self.scroll_bar.place(x=925, y=100, relheight=0.65)
        self.result_output.place(x=100, y=100)

    def update_algorithm_choice(self, choice):
        self.algorithm_choice.set(choice)

    def update_file_choice(self, choice):
        self.file_choice.set(choice)

    def clear_text(self):
        self.entry_field.delete(0, "end")
        self.result_output.delete('1.0', "end")

    """This is the main searching function
        It decides which algorithm to use depending on the user's choice
    """
    def perform_search(self):
        if self.file_choice.get() != -1:
            self.filename = PATH + self.files[self.file_choice.get()]

        self.genome_sequence = self.read_file(self.filename)
        self.query = self.entry_field.get()

        if self.genome_sequence is None:
            return
        if len(self.query) == 0:
            return None

        algorithm = self.algorithms[self.algorithm_choice.get()]
        if self.algorithm_choice.get() == 0:
            algorithm_name = "Brute Force"
        elif self.algorithm_choice.get() == 1:
            algorithm_name = "Knuth Morris Pratt"
        elif self.algorithm_choice.get() == 2:
            algorithm_name = "Rabin Karp"

        positions, count, time_elapsed, comparisons = algorithm(self.genome_sequence.upper(), self.query)
        self.result_output.delete('1.0', "end")
        text = "\n" + "Alogithm used: " + algorithm_name + "\n\n" + "File used: " + self.filename + "\n\n" + "Time elapsed: " + str(
            time_elapsed) + "\n\n" + "Number of characters in file: " + str(
            len(self.genome_sequence)) + "\n\n" + "Number of character comparisons: " + str(
            comparisons) + "\n\n" + "Count: " + str(count) + "\n\n" + "Positions: "
        if len(positions) != 0:
            positions_string = [str(element) for element in positions]

            positions_string = ", ".join(positions_string)

            string_to_append = positions_string

        else:
            string_to_append = "No occurrence"

        text += string_to_append

        self.result_output.insert('1.0', text)

    # This function tries to convert the contents of the fna file to a string
    def read_file(self, filename):
        genome = ''
        try:
            with open(filename, 'r') as f:
                for line in f:
                    if line[0] != '>':
                        genome = genome + line.rstrip()
            return genome
        except FileNotFoundError:
            messagebox.showerror("No File Selected", "Please select a file")
            return

    # This function checks if the file specified by the user exists in the Genomes Folder
    def set_file(self, file):
        if os.path.isfile(PATH + file):
            self.filename = PATH + file
            self.screen.destroy()
            messagebox.showinfo("File Found", "Using file\n" + file)
        elif os.path.isfile(file):
            self.filename = file
            self.screen.destroy()
            messagebox.showinfo("File Found", "Using file\n" + file)
        else:
            messagebox.showerror("File Not Found!",
                                 "No such file exists...\nMake sure the file is in the correct path (Genomes) and has the .fna extension")

    def on_entry_click(self, event):
        if self.entry_field.get() == self.light_text:
            self.entry_field.delete(0, "end")
            self.entry_field.insert(0, '')
            self.entry_field.config(fg='black')

    def on_focusout(self, event):
        if self.entry_field.get() == '':
            self.entry_field.insert(0, self.light_text)
            self.entry_field.config(fg='grey')

    def upload_file(self):
        self.screen = Toplevel()
        self.file_choice.set(-1)  # indicates the file is uploaded by the user
        self.screen.title("Upload File")
        label = Label(self.screen, text="Enter the filename")
        entry = Entry(self.screen, width=20, font=("Times New Roman", 16))
        enter_button = Button(self.screen, text="Enter", width=20, command=lambda: self.set_file(entry.get()))
        quit_button = Button(self.screen, text="Quit", width=20, command=self.screen.destroy)
        helper_label = Label(self.screen,
                             text="Make sure the file is in the correct path (Genomes) \nand has the .fna extension")
        label.pack(pady=10)
        entry.pack(pady=12)
        enter_button.pack(side=TOP, ipady=5, pady=5)
        quit_button.pack(side=TOP, ipady=5, pady=5)
        helper_label.pack(side=BOTTOM, pady=5)


p = ProjectApp()
