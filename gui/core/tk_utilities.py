import tkinter as tk
from tkinter import filedialog


def select_files(files_path, ext):
    """
    Open a file explorer to select files
    :param selected_files_path: list to overwrite with seletecd files
    :param label: label to update
    :param ext: valid extension
    """
    root = tk.Tk()
    root.attributes('-topmost', True)
    root.withdraw()
    file_types = [(f"{ext.upper()} files", f"*{ext}")]
    files = filedialog.askopenfilenames(parent=root, filetypes=file_types)
    root.destroy()
    add_files(files_path, files)


def select_file( ext):
    """
    Open a file explorer to select a file
    :param ext: extension file
    :return: Selected file path
    """
    root = tk.Tk()
    root.attributes('-topmost', True)
    root.withdraw()
    file_types = [(f"{ext.upper()} files", f"*{ext}")]
    file_path = filedialog.askopenfilename(filetypes=file_types)
    root.destroy()
    return file_path


def window_download_p(path):
    """
    Displays a message that the execution has been completed 
    :param path: path with the path where the protein has been downloaded
    """
    root = tk.Tk()
    root.attributes('-topmost', True)
    root.withdraw()  
    tk.messagebox.showinfo("Download", f"File saved in {path}")
    root.destroy()


def window_convert_to_mrc(path):
    """
    Displays a message that the execution has been completed 
    :param path: path where the file is saved
    """
    root = tk.Tk()
    root.attributes('-topmost', True)
    root.withdraw()  
    tk.messagebox.showinfo("Convert to mrc", f"File saved in {path}")
    root.destroy()

def window_align_to_axis(path):
    """
    Displays a message that the execution has been completed 
    :param path: path where the file is saved 
    """
    root = tk.Tk()
    root.attributes('-topmost', True)
    root.withdraw()  
    tk.messagebox.showinfo("Align to axis", f"File saved in {path}")
    root.destroy()


def window_create_helix(path):
    """
    Displays a message that the execution has been completed 
    :param path: path where the file is saved
    """
    root = tk.Tk()
    root.attributes('-topmost', True)
    root.withdraw()  
    tk.messagebox.showinfo("Create filament file", f"File saved in {path}")
    root.destroy()


def window_create_mmolecules(path):
    """
    Displays a message that the execution has been completed 
    :param path: path where the file is saved
    """
    root = tk.Tk()
    root.attributes('-topmost', True)
    root.withdraw()  
    tk.messagebox.showinfo("Create macromolecules file", f"File saved in {path}")
    root.destroy()


def window_create_membranes(path):
    """
    Displays a message that the execution has been completed 
    :param path: path where the file is saved
    """
    root = tk.Tk()
    root.attributes('-topmost', True)
    root.withdraw()  
    tk.messagebox.showinfo("Create membrane files", f"File saved in {path}")
    root.destroy()

def window_exec_app_failed():
    """
    Displays a error message that the execution has been failed
    """
    root = tk.Tk()
    root.attributes('-topmost', True)
    root.withdraw()  
    tk.messagebox.showerror("empty file list", f"Remember to add at least one file first")
    root.destroy()

def add_files(list, selected_files_path):
    """
    Adds the elements from the selected_files_path to the list.
    :param membrane_list: list containing files.
    :patam selected_files_path: list of files to be added to the list.
    """
    for file_path in selected_files_path:
        if file_path not in list:
            list.append(file_path)
    

    