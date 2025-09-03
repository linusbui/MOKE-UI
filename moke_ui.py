import os
from tkinter import *
from tkinter import ttk

import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib import rcParams
rcParams['text.usetex'] = True
rcParams['text.latex.preamble'] = r'\usepackage{{amsmath}}'

from moke_core import *
from sympy.parsing.sympy_parser import parse_expr


script_loc = os.path.dirname(os.path.abspath(__file__))
px = 1/plt.rcParams['figure.dpi']  # pixel in inches


offset_text = 0.25

PAD_SMALL = 5
PAD_MEDIUM = 10
PAD_BIG = 15


lat_k = ''
lat_g = ''
lat_h = ''

lat_xx = ''
lat_yy = ''
lat_zz = ''

lat_xy = ''
lat_yx = ''
lat_xz = ''
lat_zx = ''
lat_yz = ''
lat_zy = ''

lat_ml = ''
lat_mt = ''
lat_mlmt = ''
lat_sq = ''


# from https://blog.teclado.com/tkinter-scrollable-frames/
class ScrollableFrame(ttk.Frame):
    def __init__(self, container, *args, **kwargs):
        super().__init__(container, *args, **kwargs)
        canvas = Canvas(self, *args, **kwargs)
        scrollbar_h = ttk.Scrollbar(self, orient="horizontal", command=canvas.xview)
        self.scrollable_frame = ttk.Frame(canvas)

        self.scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(
                scrollregion=canvas.bbox("all")
            )
        )

        canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")

        canvas.configure(xscrollcommand=scrollbar_h.set)

        scrollbar_h.pack(side="bottom", fill="x")
        canvas.pack(side="left", fill="both", expand=True)


'''
def ui_routine():
    # get selections
    crystal = cryst_opt.get()
    orient = orient_opt.get()
    magn_rep = magn_reps_opt.get()
    polar = polarizatons_opt.get()


    # get latex representations based on the selections
    mo_text, latex_eps, lat_kerr, lat_eight = core_routine(crystal, orient, magn_rep, polar)

    # text for copying
    global lat_k
    lat_k = mo_text[0]
    global lat_g
    lat_g = mo_text[1]
    global lat_h
    lat_h = mo_text[2]

    global lat_xx
    lat_xx = latex_eps[0]
    global lat_yy
    lat_yy = latex_eps[1]
    global lat_zz
    lat_zz = latex_eps[2]

    global lat_xy
    lat_xy = latex_eps[3]
    global lat_yx
    lat_yx = latex_eps[4]

    global lat_xz
    lat_xz = latex_eps[5]
    global lat_zx
    lat_zx = latex_eps[6]

    global lat_yz
    lat_yz = latex_eps[7]
    global lat_zy
    lat_zy = latex_eps[8]

    global lat_ml
    lat_ml = lat_eight[0]
    global lat_mt
    lat_mt = lat_eight[1]
    global lat_mlmt
    lat_mlmt = lat_eight[2]
    global lat_sq
    lat_sq = lat_eight[3]


    # render everything
    draw_latex_MO(mo_text)
    draw_epsilon(latex_eps)
    draw_kerr(lat_kerr, lat_eight)
'''

def ui_routine():
    # get selections
    crystal = cryst_opt.get()
    orient = orient_opt.get()
    magn_rep = magn_reps_opt.get()
    polar = polarizatons_opt.get()

    #core_routine(crystal, orient, magn_rep, polar)

    # folder management
    try:
        os.makedirs(os.path.join('exports', crystal, orient, magn_rep))
        os.makedirs(os.path.join('exports', crystal, orient, magn_rep, 'p'))
        os.mkdir(os.path.join('exports', crystal, orient, magn_rep, 's'))
    except OSError:
        pass


    # read everything
    mo_show_s, mo_rot_s = get_mo(crystal, orient, magn_rep)
    mo_show = []
    mo_rot = []
    for s in mo_show_s:
        p = parse_expr(s)
        mo_show.append(p)

    for s in mo_rot_s:
        p = parse_expr(s)
        mo_rot.append(p)

    eps = get_eps(crystal, orient, mo_rot, magn_rep)
    kerr, kerr_l = get_kerr(crystal, orient, eps, magn_rep, polar)
    eight = get_eight(crystal, orient, kerr, magn_rep)

    # plotting
    draw_latex_MO(mo_show, magn_rep)
    draw_epsilon(latex_eps_comps(eps))
    draw_kerr(latex_kerr_angles(kerr), latex_eight_dir(eight))


def get_mo(crystal: str, orient: str, magn_rep: str):
    fpath_k_show = os.path.join('exports', crystal, orient, magn_rep, 'wx_k_s.txt')
    fpath_g_show = os.path.join('exports', crystal, orient, magn_rep, 'wx_g_s.txt')
    fpath_h_show = os.path.join('exports', crystal, orient, magn_rep, 'wx_h_s.txt')

    fpath_k_rot = os.path.join('exports', crystal, orient, magn_rep, 'wx_k_r.txt')
    fpath_g_rot = os.path.join('exports', crystal, orient, magn_rep,'wx_g_r.txt')
    fpath_h_rot = os.path.join('exports', crystal, orient, magn_rep, 'wx_h_r.txt')

    # all three mo tensors are created in unison, so all should be present at the same time
    if (os.path.isfile(fpath_k_show) and os.path.isfile(fpath_g_show) and os.path.isfile(fpath_h_show) and
         os.path.isfile(fpath_k_rot) and os.path.isfile(fpath_g_rot) and os.path.isfile(fpath_h_rot)):
        
        f_k_s = open(fpath_k_show, 'r')
        f_g_s = open(fpath_g_show, 'r')
        f_h_s = open(fpath_h_show, 'r')

        mo_show = [f_k_s.read(), f_g_s.read(), f_h_s.read()]

        f_k_r = open(fpath_k_rot, 'r')
        f_g_r = open(fpath_g_rot, 'r')
        f_h_r = open(fpath_h_rot, 'r')

        mo_rot = [f_k_r.read(), f_g_r.read(), f_h_r.read()]

        f_k_s.close()
        f_g_s.close()
        f_h_s.close()

        f_k_r.close()
        f_g_r.close()
        f_h_r.close()
    else:
        try:
            open(fpath_k_show, 'x')
            open(fpath_g_show, 'x')
            open(fpath_h_show, 'x')

            open(fpath_k_rot, 'x')
            open(fpath_g_rot, 'x')
            open(fpath_h_rot, 'x')
        except OSError:
            print(OSError)
            pass

        mo_show, mo_rot = core_mo(crystal, orient, magn_rep)

        f_k_s = open(fpath_k_show, 'w')
        f_g_s = open(fpath_g_show, 'w')
        f_h_s = open(fpath_h_show, 'w')

        f_k_r = open(fpath_k_rot, 'w')
        f_g_r = open(fpath_g_rot, 'w')
        f_h_r = open(fpath_h_rot, 'w')

        f_k_s.write(mo_show[0])
        f_g_s.write(mo_show[1])
        f_h_s.write(mo_show[2])

        f_k_r.write(mo_rot[0])
        f_g_r.write(mo_rot[1])
        f_h_r.write(mo_rot[2])

        f_k_s.close()
        f_g_s.close()
        f_h_s.close()

        f_k_r.close()
        f_g_r.close()
        f_h_r.close()
    
    return mo_show, mo_rot


def get_eps(crystal: str, orient: str, mo: list, magn_rep: str):
    fp_eps = os.path.join('exports', crystal, orient, magn_rep, 'w_eps.txt')

    if os.path.isfile(fp_eps):
        f_eps = open(fp_eps, 'r')
        eps_s = f_eps.read()
        eps = parse_expr(eps_s)
        f_eps.close()
    else:
        try:
            open(fp_eps, 'x')
        except OSError:
            print('error')
            pass

        eps = core_eps(mo, magn_rep)
        eps_s = str(eps)
        f_eps = open(fp_eps, 'w')
        f_eps.write(eps_s)

    return eps


def get_kerr(crystal: str, orient: str, full_eps: Matrix, magn_rep: str, polar: str):
    fp_kerr = os.path.join('exports', crystal, orient, magn_rep, polar, 'w_kerr.txt')

    if os.path.isfile(fp_kerr):
        f_kerr = open(fp_kerr, 'r')

        kerr_s = f_kerr.read()
        kerr = parse_expr(kerr_s)
        f_kerr.close()
    else:
        try:
            f_kerr = open(fp_kerr, 'x')
        except OSError:
            print('error on file creation')
            pass

        kerr = core_kerr(full_eps, magn_rep, polar)
        kerr_s = str(kerr)

        f_kerr = open(fp_kerr, 'w')
        f_kerr.write(kerr_s)
        f_kerr.close()

    return kerr, latex(kerr)


def get_eight(crystal: str, orient: str, kerr: Symbol, magn_rep: str):
    fp_eight_ml = os.path.join('exports', crystal, orient, magn_rep, 'w_ml.txt')
    fp_eight_mt = os.path.join('exports', crystal, orient, magn_rep, 'w_mt.txt')
    fp_eight_mtml = os.path.join('exports', crystal, orient, magn_rep, 'w_mtml.txt')
    fp_eight_mt2ml2 = os.path.join('exports', crystal, orient, magn_rep, 'w_mt2ml2.txt')

    if (os.path.isfile(fp_eight_ml) and os.path.isfile(fp_eight_mt) and os.path.isfile(fp_eight_mtml) and
        os.path.isfile(fp_eight_ml) and fp_eight_mt2ml2):

        f_eight_ml = open(fp_eight_ml, 'r')
        f_eight_mt = open(fp_eight_mt, 'r')
        f_eight_mtml = open(fp_eight_mtml, 'r')
        f_eight_mt2ml2 = open(fp_eight_mt2ml2, 'r')

        eight_s = [f_eight_ml.read(), f_eight_mt.read(), f_eight_mtml.read(), f_eight_mt2ml2.read()]

        eight = []
        for e in eight_s:
            eight.append(parse_expr(e))

        f_eight_ml.close()
        f_eight_mt.close()
        f_eight_mtml.close()
        f_eight_mt2ml2.close()
    else:
        try:
            f_eight_ml = open(fp_eight_ml, 'x')
            f_eight_mt = open(fp_eight_mt, 'x')
            f_eight_mtml = open(fp_eight_mtml, 'x')
            f_eight_mt2ml2 = open(fp_eight_mt2ml2, 'x')
        except OSError:
            print('error on file creation')
            pass

        eight = core_eight(kerr, magn_rep)

        f_eight_ml = open(fp_eight_ml, 'w')
        f_eight_mt = open(fp_eight_mt, 'w')
        f_eight_mtml = open(fp_eight_mtml, 'w')
        f_eight_mt2ml2 = open(fp_eight_mt2ml2, 'w')

        f_eight_ml.write(str(eight[0]))
        f_eight_mt.write(str(eight[1]))
        f_eight_mtml.write(str(eight[2]))
        f_eight_mt2ml2.write(str(eight[3]))

        f_eight_ml.close()
        f_eight_mt.close()
        f_eight_mtml.close()
        f_eight_mt2ml2.close()

    return eight


def draw_latex_MO(mo: list, magn_rep):
    M = magnetization_rep(magn_rep)
    draw_latex(wx_k, canvas_k, 0.05, 0.46, latex(mo[0]) + latex_magn_k(M))
    draw_latex(wx_g, canvas_g, 0.05, 0.46, latex(mo[1]) + latex_magn_g(M))
    draw_latex(wx_h, canvas_h, 0.009, 0.48, latex(mo[2]) + latex_magn_h(M))

    # text for copying
    global lat_k
    lat_k = latex(mo[0]) + latex_magn_k(M)
    global lat_g
    lat_g = latex(mo[1]) + latex_magn_g(M)
    global lat_h
    lat_h = latex(mo[2]) + latex_magn_h(M)

    
def draw_epsilon(latex_eps: list):
    # diagonals
    draw_latex(w_xx, canvas_xx, 0, offset_text, latex_eps[0])
    draw_latex(w_yy, canvas_yy, 0, offset_text, latex_eps[1])
    draw_latex(w_zz, canvas_zz, 0, offset_text, latex_eps[2])

    # off diags
    draw_latex(w_xy, canvas_xy, 0, offset_text, latex_eps[3])
    draw_latex(w_yx, canvas_yx, 0, offset_text, latex_eps[4])

    draw_latex(w_xz, canvas_xz, 0, offset_text, latex_eps[5])
    draw_latex(w_zx, canvas_zx, 0, offset_text, latex_eps[6])

    draw_latex(w_yz, canvas_yz, 0, offset_text, latex_eps[7])
    draw_latex(w_zy, canvas_zy, 0, offset_text, latex_eps[8])

    global lat_xx
    lat_xx = latex_eps[0]
    global lat_yy
    lat_yy = latex_eps[1]
    global lat_zz
    lat_zz = latex_eps[2]

    global lat_xy
    lat_xy = latex_eps[3]
    global lat_yx
    lat_yx = latex_eps[4]

    global lat_xz
    lat_xz = latex_eps[5]
    global lat_zx
    lat_zx = latex_eps[6]

    global lat_yz
    lat_yz = latex_eps[7]
    global lat_zy
    lat_zy = latex_eps[8]


def draw_kerr(kerr: str, eightdir: list):
    # kerr angles
    draw_latex(w_kerr, canvas_kerr, 0, offset_text, kerr)

    # eightdir
    draw_latex(w_ml, canvas_ml, 0, offset_text, eightdir[0])
    draw_latex(w_mt, canvas_mt, 0, offset_text, eightdir[1])
    draw_latex(w_mtml, canvas_mtml, 0, offset_text, eightdir[2])
    draw_latex(w_mt2ml2, canvas_mt2ml2, 0, offset_text, eightdir[3])

    global lat_ml
    lat_ml = eightdir[0]
    global lat_mt
    lat_mt = eightdir[1]
    global lat_mlmt
    lat_mlmt = eightdir[2]
    global lat_sq
    lat_sq = eightdir[3]


# subplot wx
def draw_latex(wx, canv, x, y, latstr):
    wx.clear()
    wx.axis('off')
    wx.text(x, y, '$' + latstr + '$', fontsize=12)
    canv.draw()


# subplot wx
# cached draw_latex, 
def draw_latex_c(wx, canv, x, y, latstr, fname_):
    wx.clear()
    wx.axis('off')
    wx.text(x, y, '$' + latstr + '$', fontsize=12)
    canv.draw()

    fname = fname_ + '.txt'
    fpath = os.path.join('exports', fname)
    if not os.path.isfile(fpath):
        open(fpath, 'x')
    f = open(fpath, 'w')
    f.write(latstr)


def magn_reps_observer(*args):
    magn_rep = magn_reps_opt.get()
    M = magnetization_rep(magn_rep)
    text = 'M = ' + latex(M)
    draw_latex(wx_magn, canvas_magn, 0.1, 0.45, text)

    M = magnetization_type('polar')
    text_type = '= ' + latex(M)

    draw_latex(wx_type, canvas_type, 0, 0.45, text_type)



'''
def magn_types_observer(*args):
    magn_types = magn_types_opt.get()
    M = magnetization_type(magn_types)
    text = '= ' + latex(M)
    draw_latex(wx_type, canvas_type, 0, 0.45, text)
'''


def setup_plot(height, width, container):
    fig = plt.Figure(figsize=(height*px, width*px), dpi=100, frameon=False)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01) 
    wx = fig.add_subplot(111)
    wx.axis('off')
    canvas = FigureCanvasTkAgg(fig, master=container)
    canvas.get_tk_widget().pack()
    canvas._tkcanvas.pack()
    
    return wx, canvas


root = Tk()
root.title('MOKE UI')
#root.config(bg='skyblue')


# copy to clipboard
def to_clipboard(text: str):
    root.clipboard_clear()
    root.clipboard_append('$' + text + '$')
    root.update()


# settings widget
settings_frame = Frame(root, width=300, height=500)
settings_frame.pack(fill = 'x', expand = True, padx = PAD_SMALL, pady = PAD_SMALL, side=LEFT)

# crystal settings
cryst_sett = LabelFrame(settings_frame, text = 'Crystal Settings')
cryst_sett.pack(fill = 'x', expand = True, padx = PAD_MEDIUM, pady = PAD_MEDIUM)
Label(cryst_sett, text = 'Crystal structure').grid(row = 1, column = 0, sticky = E, pady = 2)
Label(cryst_sett, text = 'Surface orientation').grid(row = 2, column = 0, sticky = E, pady = 2)

crystals = ['monoclininc', 'orthorhombic', 'tetragonal', 'tetragonal2', 'trigonal', 'trigonal2', 'hexagonal', 'hexagonal2',
            'cubic', 'cubic2']
cryst_opt = StringVar(value='	', master = cryst_sett)
OptionMenu(cryst_sett, cryst_opt, *crystals).grid(row = 1, column = 1, sticky = E, pady = 2, padx = (0,PAD_SMALL))

orientations = ['001', '011', '111']
orient_opt = StringVar(value='	', master = cryst_sett)
OptionMenu(cryst_sett, orient_opt, *orientations).grid(row = 2, column = 1, sticky = E, pady = 2, padx = (0,PAD_SMALL))


# magnetization settings
magn_sett = LabelFrame(settings_frame, text = 'Magnetization Settings')
magn_sett.pack(fill = 'x', expand = True, padx = PAD_MEDIUM, pady = PAD_MEDIUM)
Label(magn_sett, text = 'Representation').grid(row = 1, column = 0, sticky = E, pady = 2)
#clab2 = Label(magn_sett, text = 'Orientation').grid(row = 2, column = 0, sticky = E, pady = 2)

magn_reps = ['letters', 'numbers']
magn_reps_opt = StringVar(value='	')
magn_reps_opt.trace_add('write', magn_reps_observer)
OptionMenu(magn_sett, magn_reps_opt, *magn_reps).grid(row = 1, column = 1, sticky = E, pady = 2, padx = (0,PAD_SMALL))


# magnetization display
#magn_container = Frame(settings_frame, width = 400, height = 200, bg = 'grey').pack(padx = PAD_MEDIUM, pady = PAD_MEDIUM, fill = 'x')
magn_container1 = Frame(magn_sett, width = 100, height = 100, bg = 'white')
magn_container1.grid(row = 3, column = 0, sticky = E, pady = PAD_SMALL, padx = (PAD_SMALL, 0))
magn_container2 = Frame(magn_sett, width = 100, height = 100, bg = 'white')
magn_container2.grid(row = 3, column = 1, sticky = E, pady = PAD_SMALL, padx = (0,PAD_SMALL))

wx_magn, canvas_magn = setup_plot(100, 100, magn_container1)
wx_type, canvas_type = setup_plot(100, 100, magn_container2)


# polarization settings
other_sett = LabelFrame(settings_frame, text = 'Other')
other_sett.pack(fill='x', expand=True, padx = PAD_MEDIUM, pady = PAD_MEDIUM)
Label(other_sett, text='Polarization').grid(row=0, column=0, sticky=E)

polarizatons = ['s', 'p']
polarizatons_opt = StringVar(value='', master=other_sett)
OptionMenu(other_sett, polarizatons_opt, *polarizatons).grid(row = 0, column = 1, sticky = E, pady = 2, padx = (0, PAD_SMALL))


Button(settings_frame, text='Calculate', command=ui_routine).pack()


# tabs for displaying various things
tabs_container = Frame()
tabs_container.pack(side=RIGHT, padx=PAD_SMALL, pady=PAD_SMALL)

tabs = ttk.Notebook(tabs_container)
tab1 = ttk.Frame(tabs)
tab2 = ttk.Frame(tabs)
tab3 = ttk.Frame(tabs)
tab4 = ttk.Frame(tabs)
tab5 = ttk.Frame(tabs)
tabs.add(tab1, text='MO Tensors')
tabs.add(tab2, text='Eps Tensor')
tabs.add(tab3, text='Kerr Angles')
tabs.add(tab4, text='Plotting')
tabs.add(tab5, text='Simplification')
tabs.pack()


# MO Tab
k_container = LabelFrame(tab1, text='Basic Tensor K')
k_container.grid(row=0, column=0, padx=(PAD_MEDIUM, 0), pady=PAD_MEDIUM)
Button(tab1, text='cpy', command=lambda: to_clipboard(lat_k)).grid(row=0, column=1, padx=PAD_SMALL)
wx_k, canvas_k = setup_plot(250, 150, k_container)


g_container = LabelFrame(tab1, text='Basic Tensor G')
g_container.grid(row=0, column=2, pady=PAD_MEDIUM)
Button(tab1, text='cpy', command=lambda: to_clipboard(lat_g)).grid(row=0, column=4, padx=PAD_SMALL, pady=PAD_SMALL)
wx_g, canvas_g = setup_plot(330, 150, g_container)


h_container = LabelFrame(tab1, text='Basic Tensor H')
h_container.grid(row=1, columnspan=3, padx=(PAD_MEDIUM, 0), pady=PAD_MEDIUM)
Button(tab1, text='cpy', command=lambda: to_clipboard(lat_h)).grid(row=1, column=4, padx=PAD_SMALL, pady=PAD_SMALL)
wx_h, canvas_h = setup_plot(630, 270, h_container)


# Eps tab
# diags
diag_container = LabelFrame(tab2, text='Diagonal elements')
diag_container.grid(row=0, column=0, sticky=W)

diag_button_cont = Frame(tab2)
diag_button_cont.grid(row=0, column=1)
Button(diag_button_cont, text='cpy', command=lambda: to_clipboard(lat_xx)).grid(row=0, padx=PAD_SMALL, pady=PAD_SMALL)
Button(diag_button_cont, text='cpy', command=lambda: to_clipboard(lat_yy)).grid(row=1, padx=PAD_SMALL, pady=PAD_SMALL)
Button(diag_button_cont, text='cpy', command=lambda: to_clipboard(lat_zz)).grid(row=2, padx=PAD_SMALL, pady=PAD_SMALL)

scroll_container_diag = ScrollableFrame(diag_container, height = 150, width = 750)
scroll_container_diag.grid(row = 0, column = 0)


exx_container = Frame(scroll_container_diag.scrollable_frame)
exx_container.grid(row=0)
eyy_container = Frame(scroll_container_diag.scrollable_frame)
eyy_container.grid(row=1)
ezz_container = Frame(scroll_container_diag.scrollable_frame)
ezz_container.grid(row=2)

w_xx, canvas_xx = setup_plot(750, 50, exx_container)
w_yy, canvas_yy = setup_plot(750, 50, eyy_container)
w_zz, canvas_zz = setup_plot(750, 50, ezz_container)


# off diags
off_container = LabelFrame(tab2, text='Off-diagonal elements')
off_container.grid(row=1, column=0, sticky=W)

off_button_cont = Frame(tab2)
off_button_cont.grid(row=1, column=1)
Button(off_button_cont, text='cpy', command=lambda: to_clipboard(lat_xy)).grid(row=0, padx=PAD_SMALL, pady=PAD_SMALL)
Button(off_button_cont, text='cpy', command=lambda: to_clipboard(lat_yx)).grid(row=1, padx=PAD_SMALL, pady=PAD_SMALL)
Button(off_button_cont, text='cpy', command=lambda: to_clipboard(lat_xz)).grid(row=2, padx=PAD_SMALL, pady=PAD_SMALL)
Button(off_button_cont, text='cpy', command=lambda: to_clipboard(lat_zx)).grid(row=3, padx=PAD_SMALL, pady=PAD_SMALL)
Button(off_button_cont, text='cpy', command=lambda: to_clipboard(lat_yz)).grid(row=4, padx=PAD_SMALL, pady=PAD_SMALL)
Button(off_button_cont, text='cpy', command=lambda: to_clipboard(lat_zy)).grid(row=5, padx=PAD_SMALL, pady=PAD_SMALL)


scroll_container_off = ScrollableFrame(off_container, width = 750, height=300)
scroll_container_off.pack(expand=True)

exy_container = Frame(scroll_container_off.scrollable_frame)
exy_container.grid(row=3)
eyx_container = Frame(scroll_container_off.scrollable_frame)
eyx_container.grid(row=4)

w_xy, canvas_xy = setup_plot(1600, 50, exy_container)
w_yx, canvas_yx = setup_plot(1600, 50, eyx_container)


exz_container = Frame(scroll_container_off.scrollable_frame)
exz_container.grid(row=PAD_SMALL)
ezx_container = Frame(scroll_container_off.scrollable_frame)
ezx_container.grid(row=6)

w_xz, canvas_xz = setup_plot(1600, 50, exz_container)
w_zx, canvas_zx = setup_plot(1600, 50, ezx_container)


eyz_container = Frame(scroll_container_off.scrollable_frame)
eyz_container.grid(row=7)
ezy_container = Frame(scroll_container_off.scrollable_frame)
ezy_container.grid(row=8)

w_yz, canvas_yz = setup_plot(1600, 50, eyz_container)
w_zy, canvas_zy = setup_plot(1600, 50, ezy_container)


# kerr angles
full_kerr_container = LabelFrame(tab3, text='Full Kerr angles')
full_kerr_container.grid(row=0, column=0)
Button(tab3, text='cpy', command=lambda: to_clipboard(kerr_angles)).grid(row=0, column=1, padx=PAD_SMALL, pady=PAD_SMALL)

scroll_container_kerr = ScrollableFrame(full_kerr_container, height=50, width = 750)
scroll_container_kerr.pack(expand=True)

kerr_container = Frame(scroll_container_kerr.scrollable_frame)
kerr_container.grid(row=0)

w_kerr, canvas_kerr = setup_plot(1600, 50, kerr_container)


# eight directional methods
eight_dir_container = LabelFrame(tab3, text='MOKE Contributions', height=100, width=100)
eight_dir_container.grid(row=1, column=0)

eight_dir_buttons = Frame(tab3)
eight_dir_buttons.grid(row=1, column=1)

Button(eight_dir_buttons, text='cpy', command=lambda: to_clipboard(lat_ml)).grid(row = 0, padx=PAD_SMALL, pady=PAD_SMALL)
Button(eight_dir_buttons, text='cpy', command=lambda: to_clipboard(lat_mt)).grid(row = 1, padx=PAD_SMALL, pady=PAD_SMALL)
Button(eight_dir_buttons, text='cpy', command=lambda: to_clipboard(lat_mlmt)).grid(row = 2, padx=PAD_SMALL, pady=PAD_SMALL)
Button(eight_dir_buttons, text='cpy', command=lambda: to_clipboard(lat_sq)).grid(row = 3, padx=PAD_SMALL, pady=PAD_SMALL)

scroll_container_eight = ScrollableFrame(eight_dir_container, width = 750, height = 200)
scroll_container_eight.pack(expand=True)

ml_container = Frame(scroll_container_eight.scrollable_frame)
ml_container.grid(row=2)
w_ml, canvas_ml = setup_plot(1600, 50, ml_container)

mt_container = Frame(scroll_container_eight.scrollable_frame)
mt_container.grid(row=3)
w_mt, canvas_mt = setup_plot(1600, 50, mt_container)

mtml_container = Frame(scroll_container_eight.scrollable_frame)
mtml_container.grid(row=4)
w_mtml, canvas_mtml = setup_plot(1600, 50, mtml_container)

mt2ml2_container = Frame(scroll_container_eight.scrollable_frame)
mt2ml2_container.grid(row=PAD_SMALL)
w_mt2ml2, canvas_mt2ml2 = setup_plot(1600, 50, mt2ml2_container)


# plotting
plot_ellip = LabelFrame(tab4, text='Kerr Ellipticity', height=100, width=100)
plot_rot = LabelFrame(tab4, text='Kerr Rotation', height=100, width=100)
plot_ellip.pack(side=LEFT)
plot_rot.pack(side=RIGHT)


try:
    os.makedirs('exports')
except OSError:
    pass

root.mainloop()