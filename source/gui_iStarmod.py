
import sys
import os
import re
from dataclasses import dataclass
from typing import List, Optional, Tuple
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import threading

# Matplotlib for previews
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import matplotlib.pyplot as plt 

from pathlib import Path
from queue import Queue, Empty

# -------- .sm round-trip friendly parser ----------

@dataclass
class Line:
    raw: str

@dataclass
class BlankLine(Line):
    pass

@dataclass
class CommentLine(Line):
    pass

@dataclass
class KeyValueLine(Line):
    key: str
    value: str
    inline_comment: str  # includes leading " #" if present

SM_LINE_RE = re.compile(r"""
    ^\s*
    (?P<key>[A-Za-z0-9_]+)           # KEY
    \s*=\s*
    (?P<value>.*?)                   # value (non-greedy, can be empty)
    (?P<comment>\s+\#.*)?            # optional inline comment (starts with #)
    \s*$
""", re.VERBOSE)

def parse_sm(text: str) -> List[Line]:
    lines: List[Line] = []
    for raw in text.splitlines():
        if not raw.strip():
            lines.append(BlankLine(raw=raw))
            continue
        if raw.strip().startswith("#"):
            lines.append(CommentLine(raw=raw))
            continue
        m = SM_LINE_RE.match(raw)
        if m:
            key = m.group("key")
            value = m.group("value").rstrip()
            if '#' not in value:
                comment = m.group("comment") or ""
            else:
                comment = m.group("comment") or ""
                comment = value + comment
                value = ""
            lines.append(KeyValueLine(raw=raw, key=key, value=value, inline_comment=comment))
        else:
            lines.append(CommentLine(raw=raw))
    return lines

def serialize_sm(lines: List[Line]) -> str:
    out = []
    for ln in lines:
        if isinstance(ln, KeyValueLine):
            val = ln.value + "  "
            com = ln.inline_comment or ""
            out.append(f"{ln.key} = {val}{com}")
        else:
            out.append(ln.raw)
    return "\n".join(out) + ("\n" if out and not out[-1].endswith("\n") else "")

def collect_keyvalue_indices(lines: List[Line]) -> List[Tuple[int, KeyValueLine]]:
    return [(i, ln) for i, ln in enumerate(lines) if isinstance(ln, KeyValueLine)]

def kv_dict(lines: List[Line]) -> dict:
    d = {}
    for ln in lines:
        if isinstance(ln, KeyValueLine):
            d[ln.key] = ln.value
    return d

def kv_multi(lines: List[Line], key: str) -> List[str]:
    return [ln.value for ln in lines if isinstance(ln, KeyValueLine) and ln.key == key]

# -------- GUI ----------
from iStarmod_alg_class import iStarmodGUI

class SMEditorApp(tk.Tk, iStarmodGUI):
    def __init__(self, path: Optional[str] = None):
        super().__init__()
        self.title("iSTARMOD GUI")
        self.geometry("1400x800")
        self.minsize(1000, 600)
        
        s=ttk.Style()
        print (s.theme_names())

        self.file_path: Optional[str] = None
        self.lines: List[Line] = []
        self.kv_indices: List[Tuple[int, KeyValueLine]] = []

        self._build_menu()
        self._build_panes()

        if path:
            self.open_file(path)
        self.worker = None
        self.plot_queue = Queue(10)

    def _build_menu(self):
        menubar = tk.Menu(self)
        filemenu = tk.Menu(menubar, tearoff=0)
        filemenu.add_command(label="Open…", command=self.menu_open, accelerator="Ctrl+O")
        filemenu.add_command(label="Save", command=self.menu_save, accelerator="Ctrl+S")
        filemenu.add_command(label="Save As…", command=self.menu_save_as, accelerator="Ctrl+Shift+S")
        filemenu.add_separator()
        filemenu.add_command(label="Exit", command=self.menu_exit)
        menubar.add_cascade(label="File", menu=filemenu)

        viewmenu = tk.Menu(menubar, tearoff=0)
        viewmenu.add_command(label="preview PIX_EXCL", command=self.run_preview)
        viewmenu.add_command(label="Clear", command=self.clear_preview)
        menubar.add_cascade(label="Preview", menu=viewmenu)
        runmenu = tk.Menu(menubar, tearoff = 0)
        runmenu.add_command(label="iSTARMOD", command=self.menu_istarmod, accelerator="Ctrl+R")
        runmenu.add_command(label="Test View", command=self.menu_unittestreadFITS, accelerator="Ctrl+Shift+E")
        # runmenu.add_command(label="Test starot", command=self.menu_unitteststarot, accelerator="Ctrl+T")
        menubar.add_cascade(label="Run", menu=runmenu)

        self.config(menu=menubar)
        # Shortcuts
        self.bind_all("<Control-o>", lambda e: self.menu_open())
        self.bind_all("<Control-s>", lambda e: self.menu_save())
        self.bind_all("<Control-S>", lambda e: self.menu_save_as())
        self.bind_all("<Control-r>", lambda e: self.menu_istarmod())
        self.bind_all("<Control-E>", lambda e: self.menu_unittestreadFITS())
        # self.bind_all("<Control-t>", lambda e: self.menu_unitteststarot())
        

    def _build_panes(self):
        paned = ttk.PanedWindow(self, orient=tk.HORIZONTAL)
        paned.pack(fill=tk.BOTH, expand=True)

        left = ttk.Frame(paned)
        paned.add(left, weight=2)

        self.right = ttk.Frame(paned)
        paned.add(self.right, weight=8)

        self._build_leftpane(left)

        self._build_rightpane()
        
        self.stop_button.config(state="disabled")


    def _build_leftpane(self,left):
        # Left content
        top = ttk.Frame(left)
        top.pack(side=tk.TOP, fill=tk.X)
        self.path_label = ttk.Label(top, text="No file loaded", anchor="w")
        self.path_label.pack(side=tk.LEFT, padx=4, pady=6)

        container = ttk.Frame(left)
        container.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=4, pady=(0,4))

        # columns = ("key", "value", "comment")
        columns = ("key", "value")
        self.tree = ttk.Treeview(container, columns=columns, show="headings", selectmode="browse")
        self.tree.heading("key", text="Key")
        self.tree.heading("value", text="Value")
        # self.tree.heading("comment", text="Inline Comment")
        self.tree.column("key", width=150, anchor="w")
        self.tree.column("value", width=400, anchor="w")
        # self.tree.column("comment", width=280, anchor="w")

        vsb = ttk.Scrollbar(container, orient="vertical", command=self.tree.yview)
        hsb = ttk.Scrollbar(container, orient="horizontal", command=self.tree.xview)
        self.tree.configure(yscroll=vsb.set, xscroll=hsb.set)
        self.tree.grid(row=0, column=0, sticky="nsew")
        vsb.grid(row=0, column=1, sticky="ns")
        hsb.grid(row=1, column=0, sticky="ew")
        container.rowconfigure(0, weight=1)
        container.columnconfigure(0, weight=0)

        editor = ttk.LabelFrame(left, text="Edit Selected Entry")
        editor.pack(side=tk.BOTTOM, fill=tk.X, padx=2, pady=(0,4))

        self.key_var = tk.StringVar()
        self.val_var = tk.StringVar()
        self.com_var = tk.StringVar()

        ttk.Label(editor, text="Key:").grid(row=0, column=0, padx=2, pady=6, sticky="w")
        self.key_entry = ttk.Entry(editor, textvariable=self.key_var, width=20)
        self.key_entry.grid(row=0, column=0, padx=4, pady=6, sticky="e")

        ttk.Label(editor, text="Value:").grid(row=1, column=0, padx=2, pady=4, sticky="w")
        self.val_entry = ttk.Entry(editor, textvariable=self.val_var, width=50)
        self.val_entry.grid(row=1, column=0, padx=2, pady=4, sticky="e")

        ttk.Label(editor, text="Inline comment:").grid(row=2, column=0, padx=2, pady=6, sticky="w")
        self.com_entry = ttk.Entry(editor, textvariable=self.com_var, width=75)
        self.com_entry.grid(row=3, column=0, padx=2, pady=6, sticky="w")

        btnframe = ttk.Frame(editor)
        btnframe.grid(row=4, column=0, padx=2, pady=6, sticky="w")
        ttk.Button(btnframe, text="Apply", command=self.apply_changes).pack(side=tk.LEFT, padx=2)
        ttk.Button(btnframe, text="Add New", command=self.add_new).pack(side=tk.LEFT, padx=2)
        ttk.Button(btnframe, text="Delete", command=self.delete_selected).pack(side=tk.LEFT, padx=2)
        ttk.Button(btnframe, text="SaveChanges", command=self.menu_save).pack(side=tk.LEFT, padx=2)
        
        # fake_button = ttk.Button(btnframe, text= "", command= None)
        # fake_button.pack(side=tk.LEFT, fill= "none", padx=2)
        # # fake_button.pack
        # # fake_button.config(state=tk.DISABLED)
        # fake_button.config(state=tk.HIDDEN)

        self.browse_button = ttk.Button(btnframe, text= "Browse...", command= self.on_browse_path)
        self.browse_button.pack(side=tk.RIGHT, padx=2)
        self.browse_button.config(state=tk.DISABLED)

        self.tree.bind("<<TreeviewSelect>>", self.on_select)
        

    def _build_rightpane(self):
        preview_controls = ttk.Frame(self.right)
        preview_controls.pack(side=tk.TOP, fill=tk.X, padx=4, pady=(8,4))

        ttk.Button(preview_controls, text="Preview PIX_EXCL", command=self.run_preview).pack(side=tk.LEFT, padx=4)
        ttk.Button(preview_controls, text="Clear", command=self.clear_preview).pack(side=tk.LEFT, padx=4)
        self.stop_button = ttk.Button(preview_controls, text="Stop iStarmod exec", command=self.stop_exec)
        self.stop_button.pack(side=tk.LEFT, padx=4)

        self.preview_info = ttk.Label(preview_controls, text="Preview panel: figures computed from current parameters")
        self.preview_info.pack(side=tk.LEFT, padx=12)

        self.figframe = ttk.Frame(self.right)
        self.figframe.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=8, pady=8)

        self.figure = Figure(figsize=(5, 4), dpi=100)
        self.ax = self.figure.add_subplot(111)
        self.ax.set_title("No preview yet")
        self.ax.set_xlabel("Wavelength (arbitrary units)")
        self.ax.set_ylabel("Flux (arbitrary units)")
        self.canvas = FigureCanvasTkAgg(self.figure, master=self.figframe)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self.toolbar = NavigationToolbar2Tk(self.canvas, self.figframe)
        self.toolbar.update()
        self.toolbar.pack(side=tk.TOP, fill=tk.X)

        self.outcome = ttk.LabelFrame(self.right, text="Iteration Results")
        self.outcome.pack(side=tk.BOTTOM, fill=tk.X, padx=2, pady=(0,4))

        self.niter_var   = tk.StringVar()
        self.resid_var   = tk.StringVar()
        self.prm_rad_var = tk.StringVar()
        self.prm_rot_var = tk.StringVar()
        self.prm_wgt_var = tk.StringVar()
        self.sec_rad_var = tk.StringVar()
        self.sec_rot_var = tk.StringVar()
        self.sec_wgt_var = tk.StringVar()
        self.ter_rad_var = tk.StringVar()
        self.ter_rot_var = tk.StringVar()
        self.ter_wgt_var = tk.StringVar()
        self.EW1_var     = tk.StringVar()
        self.EW2_var     = tk.StringVar()
        self.EW3_var     = tk.StringVar()

        ttk.Label(self.outcome, text="N iter:").grid(row=0, column=0, padx=2, pady=6, sticky="w")
        self.niter_entry = ttk.Entry(self.outcome, textvariable=self.niter_var, width=10)
        self.niter_entry.grid(row=0, column=1, padx=2, pady=4, sticky="w")

        ttk.Label(self.outcome, text="Residuals:").grid(row=0, column=2, padx=2, pady=4, sticky="e")
        self.resid_entry = ttk.Entry(self.outcome, textvariable=self.resid_var, width=20)
        self.resid_entry.grid(row=0, column=3, padx=4, pady=6, sticky="w")

        ttk.Label(self.outcome, text="PRIMARY:  ").grid(row=1, column=1, padx=2, pady=4, sticky="w")
        ttk.Label(self.outcome, text="SECONDARY:").grid(row=1, column=2, padx=2, pady=4, sticky="w")
        ttk.Label(self.outcome, text="TERTIARY: ").grid(row=1, column=3, padx=2, pady=4, sticky="w")
        ttk.Label(self.outcome, text="Vel. Rad :").grid(row=2, column=0, padx=2, pady=4, sticky="w")
        ttk.Label(self.outcome, text="Vel. Rot :").grid(row=3, column=0, padx=2, pady=4, sticky="w")
        ttk.Label(self.outcome, text="Weight :  ").grid(row=4, column=0, padx=2, pady=4, sticky="w")
        
        self.prm_rad_entry = ttk.Entry(self.outcome, textvariable=self.prm_rad_var, width=20)
        self.prm_rad_entry.grid(row=2, column=1, padx=2, pady=6, sticky="e")
        self.sec_rad_entry = ttk.Entry(self.outcome, textvariable=self.sec_rad_var, width=20)
        self.sec_rad_entry.grid(row=2, column=2, padx=2, pady=6, sticky="e")
        self.ter_rad_entry = ttk.Entry(self.outcome, textvariable=self.ter_rad_var, width=20)
        self.ter_rad_entry.grid(row=2, column=3, padx=2, pady=6, sticky="e")
        
        self.prm_rot_entry = ttk.Entry(self.outcome, textvariable=self.prm_rot_var, width=20)
        self.prm_rot_entry.grid(row=3, column=1, padx=2, pady=6, sticky="e")
        self.sec_rot_entry = ttk.Entry(self.outcome, textvariable=self.sec_rot_var, width=20)
        self.sec_rot_entry.grid(row=3, column=2, padx=2, pady=6, sticky="e")
        self.ter_rot_entry = ttk.Entry(self.outcome, textvariable=self.ter_rot_var, width=20)
        self.ter_rot_entry.grid(row=3, column=3, padx=2, pady=6, sticky="e")

        self.prm_wgt_entry = ttk.Entry(self.outcome, textvariable=self.prm_wgt_var, width=20)
        self.prm_wgt_entry.grid(row=4, column=1, padx=2, pady=6, sticky="w")
        self.sec_wgt_entry = ttk.Entry(self.outcome, textvariable=self.sec_wgt_var, width=20)
        self.sec_wgt_entry.grid(row=4, column=2, padx=2, pady=6, sticky="w")
        self.ter_wgt_entry = ttk.Entry(self.outcome, textvariable=self.ter_wgt_var, width=20)
        self.ter_wgt_entry.grid(row=4, column=3, padx=2, pady=6, sticky="w")

        self.partial_results = []
        self.partial_results.append([self.niter_var, 0, 0])
        self.partial_results.append([self.prm_rad_var, self.prm_rot_var, self.prm_wgt_var])
        self.partial_results.append([self.sec_rad_var, self.sec_rot_var, self.sec_wgt_var])
        self.partial_results.append([self.ter_rad_var, self.ter_rot_var, self.ter_wgt_var])

        ttk.Label(self.outcome, text="Equivalent Widths (EW):").grid(row=0, column=9, padx=2, pady=4, sticky="w")
        ttk.Label(self.outcome, text="  PRIMARY:").grid(row=2, column=7, padx=2, pady=4, sticky="e")
        ttk.Label(self.outcome, text="SECONDARY:").grid(row=3, column=7, padx=2, pady=4, sticky="e")
        ttk.Label(self.outcome, text=" TERTIARY:").grid(row=4, column=7, padx=2, pady=4, sticky="e")
        self.EW1_entry = ttk.Entry(self.outcome, textvariable=self.EW1_var, width=20)
        self.EW1_entry.grid(row=2, column=9, padx=2, pady=6, sticky="e")
        self.EW2_entry = ttk.Entry(self.outcome, textvariable=self.EW2_var, width=20)
        self.EW2_entry.grid(row=3, column=9, padx=2, pady=6, sticky="e")
        self.EW3_entry = ttk.Entry(self.outcome, textvariable=self.EW3_var, width=20)
        self.EW3_entry.grid(row=4, column=9, padx=2, pady=6, sticky="e")
        
        return

    # ----- File ops -----
    def menu_open(self):
        path = filedialog.askopenfilename(
            title="Open .sm file",
            filetypes=[("StarMod files", "*.sm"), ("All files", "*.*")],
        )
        if path:
            self.open_file(path)

    def on_browse_path(self):
        key = self.key_var.get()
        value = self.val_var.get()
        # key = getattr(self, "selected_key", None)
        if not key:
            return
        
        current = self.val_var.get().strip()

        # Heuristic: directory vs file dialog
        if "_PATH" in key:
            chosen = filedialog.askdirectory(initialdir=current or None, title=f"Select folder for {key}")
        elif "_NAME" in key:
            full_path = filedialog.askopenfilename(initialdir=current or None, title=f"Select file for {key}")
            chosen = Path(full_path).name
            print(chosen)

        if not chosen:
            return  # user cancelled

        ########### Update UI + store  #################
        ## Following the CONOPS, once set the value in the edit box ###
        ## we have to 'Apply Changes' [Out of this function, in the GUI] ######
        self.val_var.set(chosen)
        ################################################
        return 

    def menu_exit(self):
        self._detach_canvas()
        if self.worker != None:
            self.stop_istarmod_event.set()
            self.worker.join()
            self.plot_queue.join()
        self.destroy()

    def menu_save(self):
        if not self.file_path:
            return self.menu_save_as()
        self.update_lines_from_tree()
        try:
            text = serialize_sm(self.lines)
            with open(self.file_path, "w", encoding="utf-8") as f:
                f.write(text)
            messagebox.showinfo("Saved", f"Saved to {self.file_path}")
        except Exception as e:
            messagebox.showerror("Error", f"Could not save:\n{e}")

    def menu_save_as(self):
        path = filedialog.asksaveasfilename(
            title="Save .sm file as…",
            defaultextension=".sm",
            filetypes=[("StarMod files", "*.sm"), ("All files", "*.*")],
        )
        if not path:
            return
        self.file_path = path
        self.menu_save()

    def open_file(self, path: str):
        try:
            with open(path, "r", encoding="utf-8") as f:
                text = f.read()
        except Exception as e:
            messagebox.showerror("Error", f"Could not open file:\n{e}")
            return
        self.file_path = path
        self.path_label.config(text=path)
        self.lines = parse_sm(text)
        self.kv_indices = collect_keyvalue_indices(self.lines)
        self.refresh_tree()

    def open_file_ref(self, path_ref):
        return

    # ----- Tree / editing -----
    def refresh_tree(self):
        self.tree.delete(*self.tree.get_children())
        for idx, kv in self.kv_indices:
            com = kv.inline_comment.strip() if kv.inline_comment else ""
            self.tree.insert("", "end", iid=str(idx), values=(kv.key, kv.value, com))
        self.key_var.set("")
        self.val_var.set("")
        self.com_var.set("")

    def on_select(self, event=None):
        sel = self.tree.selection()
        if not sel:
            return
        iid = sel[0]
        idx = int(iid)
        kv = None
        for i, k in self.kv_indices:
            if i == idx:
                kv = k
                break
        if kv is None:
            return
        self.key_var.set(kv.key)
        self.val_var.set(kv.value)
        self.com_var.set(kv.inline_comment.strip() if kv.inline_comment else "")
        ###############################################################################
        if "_PATH" in kv.key or "_NAME" in kv.key:
            self.browse_button.state(["!disabled"])
        else:
            self.browse_button.state([tk.DISABLED])
        ###############################################################################

    def apply_changes(self):
        sel = self.tree.selection()
        if not sel:
            messagebox.showwarning("No selection", "Select a row to modify.")
            return
        iid = sel[0]
        idx = int(iid)

        com_text = self.com_var.get().strip()
        if com_text and not com_text.startswith("#"):
            com_text = "# " + com_text

        for j, kv in self.kv_indices:
            if j == idx:
                kv.key = self.key_var.get().strip()
                kv.value = self.val_var.get()
                kv.inline_comment = (" " + com_text) if com_text else " # Empty inline comment"
                break
        self.tree.item(iid, values=(self.key_var.get().strip(), self.val_var.get(), com_text))

    def add_new(self):
        insert_pos = len(self.lines)
        sel = self.tree.selection()
        if sel:
            insert_pos = int(sel[0]) + 1

        new_kv = KeyValueLine(raw="", key="NEW_KEY", value="", inline_comment="")
        self.lines.insert(insert_pos, new_kv)
        self.kv_indices = collect_keyvalue_indices(self.lines)
        self.refresh_tree()
        for i, kv in self.kv_indices:
            if kv is new_kv:
                self.tree.selection_set(str(i))
                self.tree.see(str(i))
                self.on_select()
                break

    def delete_selected(self):
        sel = self.tree.selection()
        if not sel:
            messagebox.showwarning("No selection", "Select a row to delete.")
            return
        idx = int(sel[0])
        if messagebox.askyesno("Delete", "Remove this key-value line? This cannot be undone."):
            if 0 <= idx < len(self.lines) and isinstance(self.lines[idx], KeyValueLine):
                del self.lines[idx]
                self.kv_indices = collect_keyvalue_indices(self.lines)
                self.refresh_tree()

    def update_lines_from_tree(self):
        pass

    # ----- iSTARMOD algorithms and figure drawing-----

    def _attach_canvas(self, figure, ax1, ax2):
        """Attach a matplotlib Figure to the right-hand panel."""
        # Remember: self.figframe must be the frame where the figure lives


        self.figure = figure
        title = ax1.get_title()
        ax1.set_title(title, font= 'serif', size='8')
        label = ax1.get_xlabel()
        ax1.set_xlabel(label, font= 'serif', size='8')
        # label = ax1.get_ylabel()
        # ax1.set_ylabel(label, font= 'serif', size='8')
        title = ax2.get_title()
        ax2.set_title(title, font= 'serif', size='8')
        label = ax2.get_xlabel()
        ax2.set_xlabel(label, font= 'serif', size='8')
        label = ax2.get_ylabel()
        ax2.set_ylabel(label, font= 'serif', size='8')
        try:
            self.figure.set_constrained_layout(True)  # optional, but helps
        except Exception:
            pass
        size=(5, 3)
        self.figure.set_size_inches(*size, forward=True)
        self.canvas = FigureCanvasTkAgg(self.figure, master=self.figframe)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self.toolbar = NavigationToolbar2Tk(self.canvas, self.figframe)
        self.toolbar.update()
        self.toolbar.pack(side=tk.TOP, fill=tk.X)

        self.outcome.update()

    def _detach_canvas(self):
        """Remove current canvas/toolbar widgets cleanly."""
        if hasattr(self, "toolbar") and self.toolbar is not None:
            self.toolbar.destroy()
            self.toolbar = None
        if hasattr(self, "canvas") and self.canvas is not None:
            self.canvas.get_tk_widget().destroy()
            self.canvas = None
        self.figure = None

    def replace_figure(self, build_fn):
        """Swap the embedded figure. build_fn(fig) should create axes and return a dict of handles if needed."""
        self._detach_canvas()
        fig = Figure(figsize=(5, 4), dpi=100)       # or any layout you want
        handles = build_fn(fig)                     # create subplots/axes here
        # store any axes references you need:
        for k, v in (handles or {}).items():
            setattr(self, k, v)
        self._attach_canvas(fig)

    def menu_unittestreadFITS(self):
        
        if (len(self.figure.get_axes()) == 2):
            self.clear_preview()
                
        # self.open_file("unittestreadFITS.sm")
        import unittestreadFITS_canvas as readFITS
        self.figure, self.ax = readFITS.unittestreadFITS(self.file_path, self.figure, self.ax)
        self.canvas.draw()
        # self.canvas = FigureCanvasTkAgg(self.figure, master=self.figframe)
        # self.canvas.draw()
        # self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        # self.toolbar = NavigationToolbar2Tk(self.canvas, self.figframe)
        # self.toolbar.update()
        # self.toolbar.pack(side=tk.TOP, fill=tk.X)

        # self.outcome.update()
        return
    
    def update_partial_results_frm(self, partial_results_var):
        
        self.niter_var.set(partial_results_var[0][0])
        self.resid_var.set(partial_results_var[0][1])
        self.prm_rad_var.set(partial_results_var[1][0])
        self.prm_rot_var.set(partial_results_var[1][1])
        self.prm_wgt_var.set(partial_results_var[1][2])
        self.sec_rad_var.set(partial_results_var[2][0])
        self.sec_rot_var.set(partial_results_var[2][1])
        self.sec_wgt_var.set(partial_results_var[2][2])
        self.ter_rad_var.set(partial_results_var[3][0])
        self.ter_rot_var.set(partial_results_var[3][1])
        self.ter_wgt_var.set(partial_results_var[3][2])

        self.canvas.draw_idle()
        # cancelling_update = self.outcome.after_idle(self.outcome.update)
        self.outcome.update_idletasks()
        
    def update_final_results(self, equivalent_widths_var):

        self.EW1_var.set(equivalent_widths_var[0])
        self.EW2_var.set(equivalent_widths_var[1])
        self.EW3_var.set(equivalent_widths_var[2])

        self.canvas.draw_idle()
        self.outcome.update_idletasks()

    
    def update_canvas_from_istarmod(self, payload):
        
        f, ax1, ax2 = self._build_plot_(payload)
        self._detach_canvas()
        self._attach_canvas(f, ax1, ax2)
        self.canvas.draw_idle()
        # self.after(150)

    def poll_plot_queue(self):
        try:
            while True:
                payload = self.plot_queue.get_nowait()
                if "error" in payload:
                    print("Worker error:", payload["error"])
                else:
                    self.update_canvas_from_istarmod(payload)
        except Empty:
            pass
        self.after(50, self.poll_plot_queue)
    
    def idle(self):
        self.update_idletasks()

    def istarmod_work(self):
        self.starmod(self.file_path, self.figure, self.ax1, self.ax2, self)
        return
    
    def menu_istarmod(self):

        self.clear_preview()
        if (len(self.figure.get_axes()) == 1):
            self.figure.delaxes(self.ax)
            self.ax1 = self.figure.add_subplot(211)
            self.ax2 = self.figure.add_subplot(212, sharex = self.ax1)
            self.canvas.draw()
            # self.figure, (self.ax1,self.ax2) = plt.subplots(2, figsize = (10,7), dpi = 150, sharex = True)#, sharey = True)
        
        if(self.file_path != None and self.file_path != "unittestreadFITS.sm" ):
            self.stop_istarmod_event = threading.Event()
            self.worker = threading.Thread(target=self.istarmod_work , daemon=True)
            # self.figure, (self.ax1,self.ax2) = plt.subplots(2, figsize = (10,7), dpi = 150, sharex = True)
            self.worker.start()

        print("Here1!!!!!")
        self.canvas.draw_idle()
        # self.canvas.after()
        
        return
    
    def menu_unitteststarot(self):
        with open("unitteststarot.py") as f:
            code = f.read()
        exec(code)

    # ----- Preview logic -----
    def run_preview(self):
        if not self.lines:
            self.ax.cla()
            self.ax.set_title("No file loaded")
            self.ax.set_xlabel("Wavelength (arbitrary units)")
            self.ax.set_ylabel("Flux (arbitrary units)")
            self.canvas.draw()
            return

        d = kv_dict(self.lines)
        excl_values = kv_multi(self.lines, "PIX_EXCL")

        # Parse PIX_ZONE like: "6530 6615  wvl"
        zone = (0.0, 1.0)
        if "PIX_ZONE" in d:
            toks = d["PIX_ZONE"].split()
            if len(toks) >= 2:
                try:
                    zone = (float(toks[0]), float(toks[1]))
                except ValueError:
                    pass
        
        import numpy as np
        from iStarmod_tools import lambdaData as lambdas
        self.ax.cla()
        x = np.linspace(zone[0] - 20, zone[1] + 20, 2000)
        l = lambdas("lambdas.dat")
        for key, value in l.lambdaDataDict.items():
            flambda = float(value[0])
            if (flambda > zone[0] and flambda < zone[1]):
                mu = flambda 
        # mu = 6564.60
        # sigma = 0.6
                if(key == "HeD3"):
                    sigma = 0.45
                    y = 1.0 - 0.2 * np.exp(-0.5 * ((x - mu) / sigma) ** 2)
                elif(key == "NaD1" or key == "NaD2" or key == "CaIIH" or key == "CaIIK" or key == "Halpha"):
                    sigma = 0.9
                    y = 1.0 - 0.7 * np.exp(-0.5 * ((x - mu) / sigma) ** 2)
                else:
                    sigma = 0.7
                    y = 1.0 - 0.5 * np.exp(-0.5 * ((x - mu) / sigma) ** 2)
                self.ax.plot(x, y, linewidth=1.0, label = key)

        spans = []
        for v in excl_values:
            toks = v.split()
            if len(toks) >= 2:
                try:
                    a, b = float(toks[0]), float(toks[1])
                    spans.append((a, b))
                except ValueError:
                    pass

        # self.ax.cla()
        # self.ax.plot(x, y, linewidth=1.0)
        self.ax.set_title("Preview: Toy spectrum with zone/exclusions")
        self.ax.set_xlabel("Wavelength")
        self.ax.set_ylabel("Flux")
        
        self.ax.legend()
        self.ax.axvspan(zone[0], zone[1], alpha=0.1)
        for a, b in spans:
            self.ax.axvspan(a, b, alpha=0.15, hatch="///")

        if "N_ITER" in d:
            try:
                niter = int(d["N_ITER"].split()[0])
                self.ax.text(0.02, 0.95, f"N_ITER = {niter}", transform=self.ax.transAxes, va="top")
            except Exception:
                pass

        self.canvas.draw()

    def clear_preview(self):

        axes_nr = len(self.figure.get_axes())
        if(axes_nr == 2):
            print("delaxes")
            # self.figure.delaxes(self.figure.get_axes()[0]) #ax1
            # self.figure.delaxes(self.figure.get_axes()[0]) #ax2
            # self._build_rightpane()
            self._detach_canvas()
            self.figure = Figure(figsize=(5, 4), dpi=100)
            self.ax = self.figure.add_subplot(111)
            self.ax.set_title("No preview yet")
            self.ax.set_xlabel("Wavelength (arbitrary units)")
            self.ax.set_ylabel("Flux (arbitrary units)")
            self.canvas = FigureCanvasTkAgg(self.figure, master=self.figframe)
            self.canvas.draw()
            self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
            self.clear_results_entries()
        
        self.ax.cla()
        self.ax.set_title("Preview cleared")
        self.ax.set_xlabel("Wavelength (arbitrary units)")
        self.ax.set_ylabel("Flux (arbitrary units)")
        self.canvas.draw()
        return
    
    def clear_results_entries(self):
        self.niter_var.set("")
        self.resid_var.set("")
        self.prm_rad_var.set("")
        self.prm_rot_var.set("")
        self.prm_wgt_var.set("")
        self.sec_rad_var.set("")
        self.sec_rot_var.set("")
        self.sec_wgt_var.set("")
        self.ter_rad_var.set("")
        self.ter_rot_var.set("")
        self.ter_wgt_var.set("")
        self.EW1_var.set("")
        self.EW2_var.set("")
        self.EW3_var.set("")

    def finalize_execution(self):
        self.stop_istarmod_event.set()
        self.worker.join()
        return

    def stop_exec(self):
        self.update_idletasks()
        if self.worker.is_alive():
            if hasattr(self, "stop_istarmod_event"):
                self.finalize_execution()
                if(self.figure.get_axes()==2):
                    print("delaxes")
                    self.figure.delaxes(self.ax1)
                    self.figure.delaxes(self.ax2)
                self.clear_preview()
                print("Here2!!!")
    
    
    def _build_plot_(self, payload):

        ################ THIS FUNCTION ISOLATES PLOTTING FROM THE ALGORTIHM ##########
        ################ the payload conatains all paramters and objects needed to plot the different figures #######

        inputParams = payload["inputParams"]
        nstar       = payload["nstar"]
        _date       = payload["date"]
        _object     = payload["object"]
        lambdas     = payload["lambdas"]
        ldoObs      = payload["ldoObs"]
        ldoObs2     = payload["ldoObs2"]
        param2      = payload["param2"]
        param4      = payload["param4"]
        param5      = payload["param5"] 
        param6      = payload["param6"]
        debugging   = payload["debugging"]
        workingLambdaValues     = payload["workingLambdaValues"]
        workingDataValues13     = payload["workingDataValues3"]
        workingDataValues21     = payload["workingDataValues1"]
        workingDataValues22     = payload["workingDataValues2"]
        summ_up                 = payload["summ_up"]
        flambdavalues           = payload["flambdavalues"] 
        star1_spectrum          = payload["star1_spectrum"]
        star2_spectrum          = payload["star2_spectrum"]
        

        # These lines must/should be commented in the worker thread in order not to Starting a Matplotlib GUI outside of the main thread 
        # (And then fail)
        f, (subPlot1,subPlot2) = plt.subplots(2, figsize = (10,7), dpi = 150, sharex = True)#, sharey = True)
        # subPlot1.clear()
        # subPlot2.clear()
        if(inputParams.linespcf != 'NONE'):
            if nstar == 1:
            ## FOR DEBUGGING PURPOSES. To draw the boundaries of the zone where to search for the maximum.
                if (debugging == True):
                    subPlot1.plot(ldoObs, param2, 'b+')
                    subPlot1.axvline(x=param5, c = 'r', linestyle='dashed')
                    subPlot1.axvline(x=param6, c = 'r', linestyle='dashed')
            elif (nstar == 2):
                if (debugging == True):
                    subPlot1.axvline(x=param5, c = 'r', linestyle='dashed')
                    subPlot1.axvline(x=param6, c = 'r', linestyle='dashed')
                    subPlot1.axvspan(param5, param6, color='green', alpha=0.1)
            ###############################################################################
            ########## To draw the integration limits #####################################
            subPlot1.axvline(x=workingLambdaValues[lambdas.start], linestyle='dashed')
            subPlot1.axvline(x=workingLambdaValues[lambdas.end],   linestyle='dashed')
            if debugging == True:
                subPlot1.axvspan(workingLambdaValues[lambdas.start], 
                                    workingLambdaValues[lambdas.end]  , color='blue', alpha=0.1)
            subPlot1.axhline(y=0.0, linestyle='dashed')
            subPlot2.axvline(x=workingLambdaValues[lambdas.start], linestyle='dashed')
            subPlot2.axvline(x=workingLambdaValues[lambdas.end],   linestyle='dashed')
            ###############################################################################
            if nstar ==2:
                    subPlot1.plot(flambdavalues, star1_spectrum, 'c-', lw = 1.5)
                    subPlot1.plot(flambdavalues, star2_spectrum, 'y-', lw = 1.5)
                    if (debugging == True):## FOR DEBUGGING PURPOSES
                        if (ldoObs2 != None and param4 != None):
                            subPlot1.plot(ldoObs2, param4, 'r*') 
                        if(ldoObs != None and param2 != None):
                            subPlot1.plot(ldoObs, param2, 'b+')
                        if summ_up != None:
                            subPlot1.plot(flambdavalues, summ_up, c = 'm', linestyle = 'dashed', lw = 1.25)


        if len(workingLambdaValues) == len(workingDataValues22):
            subPlot2.plot(workingLambdaValues, workingDataValues21, 'b-', lw = 1.0)
            subPlot2.plot(workingLambdaValues, workingDataValues22, 'r-', lw = 1.0)
            subPlot1.plot(workingLambdaValues, workingDataValues13, 'g-', lw = 1.0)
        elif len(workingLambdaValues) >= len(workingDataValues22):
            subPlot2.plot(workingLambdaValues, workingDataValues21, 'b-', lw = 1.0)
            subPlot2.plot(workingLambdaValues[0:len(workingDataValues22)], workingDataValues22, 'r-', lw = 1.0)
            subPlot1.plot(workingLambdaValues, workingDataValues13, 'g-', lw = 1.0)
        elif len(workingLambdaValues) <= len(workingDataValues22):
            subPlot2.plot(workingLambdaValues, workingDataValues21, 'b-', lw = 1.0)
            subPlot2.plot(workingLambdaValues, workingDataValues22[0:len(workingLambdaValues)], 'r-', lw = 1.0)
            subPlot1.plot(workingLambdaValues, workingDataValues13, 'g-', lw = 1.0)

        f.subplots_adjust(hspace = 0.3)

        dateTitle = str('{:6.4f}'.format(_date))
        
        strTitle = str(_object) + "@" + str(dateTitle) + "_" + inputParams.linespcf 
        if (inputParams.wvl_display_range[0] != -1 and inputParams.wvl_display_range[1] != -1): 
            subPlot1.set_xlim(inputParams.wvl_display_range[0],inputParams.wvl_display_range[1])
        # # CaIIH&K
        if(inputParams.linespcf == "CaIIH" or inputParams.linespcf == "CaIIK"):
            subPlot1.set_title("Subtracted Spectrum for "+ inputParams.linespcf + r" / $v_{\rm rot}sini$ = " + str('{:7.4f}'.format(inputParams.PRMRot[0])) + " km/s", font = 'serif', size= 18)
            subPlot1.set_ylim(-1.0, inputParams.MaxFluxDisplayed_Sub)
            subPlot2.set_title(_object + "@" + dateTitle +"  &  Synthetic Spectra for " + inputParams.linespcf, font = 'serif', size= 18)
            subPlot2.set_ylim(-0.5, inputParams.MaxFluxDisplayed_Obj)
            # subPlot1.set_xlim(3945,3985)
            # subPlot2.set_xlim(3945,3985)
        # # Hbeta
        elif(inputParams.linespcf == "Hbeta"):
            subPlot1.set_title("Subtracted Spectrum for "+ r"$H\beta$" + r" / $v_{\rm rot}sini$ = " + str('{:7.4f}'.format(inputParams.PRMRot[0])) + " km/s", font = 'serif', size= 18)
            subPlot1.set_ylim(-1.0, inputParams.MaxFluxDisplayed_Sub)
            subPlot2.set_title(_object + "@" + dateTitle + "  &  Synthetic Spectra for " + r"$H\beta$", font = 'serif', size= 18)
            subPlot2.set_ylim(-0.5, inputParams.MaxFluxDisplayed_Obj)
            if (inputParams.wvl_display_range[0] == -1 and inputParams.wvl_display_range[1] == -1):
                subPlot1.set_xlim(4850,4890)
                
        # # HeD3 NaD2/NaD1
        elif(inputParams.linespcf == 'HeD3' or inputParams.linespcf == 'NaD2' or inputParams.linespcf == 'NaD1'):
            subPlot1.set_title(r"Subtracted Spectrum for "+ inputParams.linespcf + r" / $v_{\rm rot}sini$ = " + str('{:7.4f}'.format(inputParams.PRMRot[0])) + " km/s", font = 'serif', size= 18)
            subPlot1.set_ylim(-1.0, inputParams.MaxFluxDisplayed_Sub)
            subPlot2.set_title(_object + "@" + dateTitle + "  &  Synthetic Spectra for " + inputParams.linespcf, font = 'serif', size= 18)
            subPlot2.set_ylim(-0.5, inputParams.MaxFluxDisplayed_Obj)
            if (inputParams.wvl_display_range[0] == -1 and inputParams.wvl_display_range[1] == -1):
                subPlot1.set_xlim(5860,5910)

        # # Halpha
        elif(inputParams.linespcf == 'Halpha'):
            subPlot1.set_title(r"Subtracted Spectrum for "+ r"$H\alpha$" + r" / $v_{\rm rot}sini$ = " + str('{:7.4f}'.format(inputParams.PRMRot[0])) + " km/s", font = 'serif', size= 18)
            subPlot2.set_title(_object + "@" + dateTitle + "  &  Synthetic Spectra for " + r"$H\alpha$", font = 'serif', size= 18)
            subPlot2.set_ylim(0.0, inputParams.MaxFluxDisplayed_Obj)
            subPlot1.set_ylim(-1.0, inputParams.MaxFluxDisplayed_Sub)
            if (inputParams.wvl_display_range[0] == -1 and inputParams.wvl_display_range[1] == -1):
                subPlot1.set_xlim(6550,6580)

        #  CaIRT-a
        elif(inputParams.linespcf == 'CaIRT-a'):
            subPlot1.set_title(r"Subtracted Spectrum for "+ inputParams.linespcf + r" / $v_{\rm rot}sini$ = " + str('{:7.4f}'.format(inputParams.PRMRot[0])) + " km/s", font = 'serif', size= 18)
            subPlot2.set_title(_object + "@" + dateTitle +"  &  Synthetic Spectra for " + inputParams.linespcf, font = 'serif', size= 18)
            subPlot1.set_ylim(-0.75, inputParams.MaxFluxDisplayed_Sub)
            subPlot2.set_ylim(0.25,inputParams.MaxFluxDisplayed_Obj)
            if (inputParams.wvl_display_range[0] == -1 and inputParams.wvl_display_range[1] == -1):
                subPlot1.set_xlim(8465,8535)

        elif(inputParams.linespcf == 'CaIRT-b'):
            subPlot1.set_title(r"Subtracted Spectrum for "+ inputParams.linespcf + r" / $v_{\rm rot}sini$ = " + str('{:7.4f}'.format(inputParams.PRMRot[0])) + " km/s", font = 'serif', size= 18)
            subPlot2.set_title(_object + "@" + dateTitle +"  &  Synthetic Spectra for " + inputParams.linespcf, font = 'serif', size= 18)
            subPlot1.set_ylim(-0.75, inputParams.MaxFluxDisplayed_Sub)
            subPlot2.set_ylim(0.0,inputParams.MaxFluxDisplayed_Obj)
            if (inputParams.wvl_display_range[0] == -1 and inputParams.wvl_display_range[1] == -1):
                subPlot1.set_xlim(8535,8550)

        elif(inputParams.linespcf == 'CaIRT-c'):
            subPlot1.set_title(r"Subtracted Spectrum for "+ inputParams.linespcf + r" / $v_{\rm rot}sini$ = " + str('{:7.4f}'.format(inputParams.PRMRot[0])) + " km/s", font = 'serif', size= 18)
            subPlot2.set_title(_object + "@" + dateTitle +"  &  Synthetic Spectra for " + inputParams.linespcf, font = 'serif', size= 18)
            subPlot1.set_ylim(-0.5, inputParams.MaxFluxDisplayed_Sub)
            subPlot2.set_ylim(0.,inputParams.MaxFluxDisplayed_Obj)
            if (inputParams.wvl_display_range[0] == -1 and inputParams.wvl_display_range[1] == -1):
                subPlot1.set_xlim(8600,8680)

        # # #################### Paschen D
        elif(inputParams.linespcf == 'PaschenD'):
            subPlot2.set_title(_object + "@" + dateTitle +"  &  Synthetic Spectra for " + inputParams.linespcf, font = 'serif', size= 18)
            subPlot2.set_ylim(0.25, inputParams.MaxFluxDisplayed_Obj)
            subPlot1.set_title(r"Subtracted Spectrum for "+ inputParams.linespcf + r" / $v_{\rm rot}sini$ = " + str('{:7.4f}'.format(inputParams.PRMRot[0])) + " km/s", font = 'serif', size= 18)
            subPlot1.set_ylim(-0.25, inputParams.MaxFluxDisplayed_Sub)
            if (inputParams.wvl_display_range[0] == -1 and inputParams.wvl_display_range[1] == -1):
                subPlot1.set_xlim(10040,10065)

        #  ##################### HeI10830
        elif(inputParams.linespcf == 'HeI10833'):
            subPlot2.set_title(_object + "@" + dateTitle +"  &  Synthetic Spectra for " + inputParams.linespcf, font = 'serif', size= 18)
            subPlot2.set_ylim(0.25, inputParams.MaxFluxDisplayed_Obj)
            subPlot1.set_title(r"Subtracted Spectrum for "+ inputParams.linespcf + r" / $v_{\rm rot}sini$ = " + str('{:7.4f}'.format(inputParams.PRMRot[0])) + " km/s", font = 'serif', size= 18)
            subPlot1.set_ylim(-0.25, inputParams.MaxFluxDisplayed_Sub)
            if (inputParams.wvl_display_range[0] == -1 and inputParams.wvl_display_range[1] == -1):
                subPlot1.set_xlim(10820,10860)

        #  ##################### PaschenG
        elif(inputParams.linespcf == 'PaschenG'):
            subPlot2.set_title(_object + "@" + dateTitle +"  &  Synthetic Spectra for " + inputParams.linespcf, font = 'serif', size= 18)
            subPlot2.set_ylim(0.25, inputParams.MaxFluxDisplayed_Obj)
            subPlot1.set_title(r"Subtracted Spectrum for "+ inputParams.linespcf + r" / $v_{\rm rot}sini$ = " + str('{:7.4f}'.format(inputParams.PRMRot[0])) + " km/s", font = 'serif', size= 18)
            subPlot1.set_ylim(-0.25, inputParams.MaxFluxDisplayed_Sub)
            if (inputParams.wvl_display_range[0] == -1 and inputParams.wvl_display_range[1] == -1):
                subPlot1.set_xlim(10920,10955)

        #  ##################### PaschenB
        elif(inputParams.linespcf == 'PaschenB'):
            subPlot2.set_title(_object + "@" + dateTitle +"  &  Synthetic Spectra for " + inputParams.linespcf, font = 'serif', size= 18)
            subPlot2.set_ylim(0.5, inputParams.MaxFluxDisplayed_Obj)
            subPlot1.set_title(r"Subtracted Spectrum for "+ inputParams.linespcf + r" / $v_{\rm rot}sini$ = " + str('{:7.4f}'.format(inputParams.PRMRot[0])) + " km/s", font = 'serif', size= 18)
            subPlot1.set_ylim(-0.5, inputParams.MaxFluxDisplayed_Sub)
            if (inputParams.wvl_display_range[0] == -1 and inputParams.wvl_display_range[1] == -1):
                subPlot1.set_xlim(12800, 12840)
                
        else:
        # Any
            subPlot2.set_title(_object + "@" + dateTitle +"  &  Synthetic Spectra for " + inputParams.linespcf, font = 'serif', size= 18)
            subPlot2.set_ylim(0.0, inputParams.MaxFluxDisplayed_Obj)
            subPlot1.set_title(r"Subtracted Spectrum for "+ inputParams.linespcf + r" / $v_{\rm rot}sini$ = " + str('{:7.4f}'.format(inputParams.PRMRot[0])) + " km/s", font = 'serif', size= 18)
            subPlot1.set_ylim(-0.5, inputParams.MaxFluxDisplayed_Sub)
            # subPlot1.set_xlim(9326 , 9485)
            # subPlot2.set_xlim(9326 , 9485)
            
        plt.xlabel(r"$\lambda$  [${\rm \AA}$]", font = 'serif', size= 14)
        plt.ylabel("Normalised Flux", font= 'serif', size= 14)
        # plt.legend()
        # # At this point is supposed that a valid value of inputParams.resultspath has been assigned
        figure_repo = Path(inputParams.resultspath) / Path("RES") / Path(inputParams.linespcf)
        if not os.access(Path(figure_repo), os.F_OK): 
            os.mkdir(Path(figure_repo))
        figureName = os.path.join(figure_repo, strTitle)
        if nstar == 1:  
                figureName = figureName + ".png"
        elif nstar == 2:
                figureName = figureName + "_asSB2.png"
        
        figureName_p = Path(figureName)
        f.savefig(figureName_p, dpi = 300)
        
        return f, subPlot1, subPlot2


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else None
    app = SMEditorApp(path)

    def on_close():
        app.quit()   # exits mainloop cleanly
        app.destroy()

    app.protocol("WM_DELETE_WINDOW", on_close)
    app.poll_plot_queue()

    try:
        if "clam" in ttk.Style().theme_names():
            ttk.Style().theme_use("clam")
    except Exception:
        pass

    try:
        app.mainloop()
        app.plot_queue.join()
    finally:
        try:
           if app.winfo_exists():
               app.destroy()
        except Exception:
            pass
    
    return


if __name__ == "__main__":
    main()
