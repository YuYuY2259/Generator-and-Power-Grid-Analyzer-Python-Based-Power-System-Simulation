import sys
import numpy as np
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                             QHBoxLayout, QTabWidget, QLabel, QSlider, 
                             QPushButton, QFrame, QScrollArea)
from PyQt6.QtCore import Qt, QTimer
from PyQt6.QtGui import QFont

# Matplotlib integration for PyQt
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib.patches as mpatches

# --- ÖZEL TOOLBAR SINIFI ---
class CustomToolbar(NavigationToolbar):
    def __init__(self, canvas, parent, coordinates=True, on_home=None):
        super().__init__(canvas, parent, coordinates)
        self.on_home_callback = on_home

    def home(self, *args):
        super().home(*args)
        if self.on_home_callback:
            self.on_home_callback()

class GeneratorAnalyzer(QMainWindow):
    def __init__(self):
        super().__init__()
        
        self.setWindowTitle("Generator Analyzer V5.13 (Current Label Added)")
        self.setGeometry(100, 50, 1100, 850)
        self.setStyleSheet("background-color: white;")

        # --- Constants ---
        self.S_base_sys = 850e6
        self.S_gen_rated = 800e6
        self.Z_xfm_pu = complex(0.0025, 0.057) 

        # --- Tab Structure ---
        self.tabs = QTabWidget()
        self.setCentralWidget(self.tabs)
        
        self.tab_home = QWidget()
        self.tab_info = QWidget()
        self.tab_about = QWidget()
        
        self.tabs.addTab(self.tab_home, "📊 Simulation")
        self.tabs.addTab(self.tab_info, "ℹ️ Info")
        self.tabs.addTab(self.tab_about, "👤 About Us")

        # Build UI
        self.setup_home_tab()
        self.setup_info_tab()
        self.setup_about_tab()

        # Initial Simulation Update
        self.update_sim()

    # --- HOME TAB SETUP ---
    def setup_home_tab(self):
        main_layout = QVBoxLayout(self.tab_home)

        # Başlık
        header_text_label = QLabel("GENERATOR AND POWER GRID ANALYSIS")
        header_text_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        header_text_label.setFont(QFont("Arial", 24, QFont.Weight.Bold))
        header_text_label.setStyleSheet("color: #333; margin-bottom: 10px;")
        main_layout.addWidget(header_text_label)

        # İçerik Düzeni
        content_layout = QHBoxLayout()
        main_layout.addLayout(content_layout)

        # Grafik Objeleri
        self.fig = Figure(figsize=(8, 10), dpi=100)
        self.canvas = FigureCanvas(self.fig)
        self.toolbar = CustomToolbar(self.canvas, None, on_home=self.reset_simulation)
        self.toolbar.setStyleSheet("background-color: white; border: none;")

        # Sol Panel
        controls_panel = QFrame()
        controls_panel.setFixedWidth(300)
        controls_layout = QVBoxLayout(controls_panel)
        controls_layout.addSpacing(85)
        controls_layout.addWidget(self.toolbar)
        controls_layout.addSpacing(15)
        
        # Kontroller
        self.lbl_p_title = QLabel("Active Power (MW):")
        self.lbl_p_title.setFont(QFont("Arial", 10, QFont.Weight.Bold))
        self.slider_p = QSlider(Qt.Orientation.Horizontal)
        self.slider_p.setRange(0, 850); self.slider_p.setValue(750)
        self.lbl_p_val = QLabel("750 MW")
        self.slider_p.valueChanged.connect(self.update_sim)
        
        self.lbl_pf_title = QLabel("Load Power Factor:")
        self.lbl_pf_title.setFont(QFont("Arial", 10, QFont.Weight.Bold))
        self.slider_pf = QSlider(Qt.Orientation.Horizontal)
        self.slider_pf.setRange(-60, 60); self.slider_pf.setValue(int(np.degrees(np.arccos(0.9))))
        self.lbl_pf_val = QLabel("0.90 Lead"); self.lbl_pf_val.setStyleSheet("color: red;")
        self.slider_pf.valueChanged.connect(self.update_sim)

        self.lbl_xd_title = QLabel("Xg on Base (pu):")
        self.lbl_xd_title.setFont(QFont("Arial", 10, QFont.Weight.Bold))
        self.slider_xd = QSlider(Qt.Orientation.Horizontal)
        self.slider_xd.setRange(50, 250); self.slider_xd.setValue(128)
        self.lbl_xd_val = QLabel("j1.28 pu")
        self.slider_xd.valueChanged.connect(self.update_sim)

        sep_line = QLabel("---------------------------------------------------")
        sep_line.setStyleSheet("color: gray;")

        self.btn_mode = QPushButton("ADVANCED MODE: OFF")
        self.btn_mode.setCheckable(True)
        self.btn_mode.setStyleSheet("background-color: #E6E6E6; font-weight: bold; padding: 5px;")
        self.btn_mode.clicked.connect(self.toggle_mode)

        self.lbl_load_title = QLabel("Load Factor:")
        self.lbl_load_title.setFont(QFont("Arial", 10, QFont.Weight.Bold))
        self.slider_load = QSlider(Qt.Orientation.Horizontal)
        self.slider_load.setRange(50, 150); self.slider_load.setValue(100); self.slider_load.setEnabled(False)
        self.lbl_load_val = QLabel("1.0x"); self.lbl_load_val.setStyleSheet("color: gray;")
        self.slider_load.valueChanged.connect(self.update_sim)

        self.lbl_vr = QLabel("VR: -- %")
        self.lbl_vr.setStyleSheet("color: blue; font-weight: bold; font-size: 14px; margin-top: 10px;")
        
        self.lbl_loss_title = QLabel("Extra Losses:")
        self.lbl_loss_title.setFont(QFont("Arial", 10, QFont.Weight.Bold))
        self.slider_loss = QSlider(Qt.Orientation.Horizontal)
        self.slider_loss.setRange(0, 50); self.slider_loss.setValue(0); self.slider_loss.setEnabled(False)
        self.lbl_loss_val = QLabel("+0.00 pu"); self.lbl_loss_val.setStyleSheet("color: gray;")
        self.slider_loss.valueChanged.connect(self.update_sim)

        self.lbl_eff = QLabel("Eff: -- %")
        self.lbl_eff.setStyleSheet("color: green; font-weight: bold; font-size: 14px; margin-top: 5px;")

        widgets = [self.lbl_p_title, self.slider_p, self.lbl_p_val, 
                   self.lbl_pf_title, self.slider_pf, self.lbl_pf_val,
                   self.lbl_xd_title, self.slider_xd, self.lbl_xd_val,
                   sep_line, self.btn_mode,
                   self.lbl_load_title, self.slider_load, self.lbl_load_val,
                   self.lbl_vr, self.lbl_loss_title, self.slider_loss, self.lbl_loss_val,
                   self.lbl_eff]
        
        for w in widgets: controls_layout.addWidget(w)
        controls_layout.addStretch()

        # Sağ Panel
        plot_panel = QWidget()
        plot_layout = QVBoxLayout(plot_panel)
        plot_layout.addWidget(self.canvas)
        
        self.ax_phasor = self.fig.add_subplot(211) 
        self.ax_bar = self.fig.add_subplot(212)    
        self.fig.subplots_adjust(hspace=0.3)
        self.canvas.mpl_connect('scroll_event', self.on_scroll)

        content_layout.addWidget(controls_panel)
        content_layout.addWidget(plot_panel, stretch=1)

    def toggle_mode(self):
        is_on = self.btn_mode.isChecked()
        if is_on:
            self.btn_mode.setText("ADVANCED MODE: ON")
            self.btn_mode.setStyleSheet("background-color: #99FF99; font-weight: bold; padding: 5px;")
            self.slider_load.setEnabled(True); self.slider_loss.setEnabled(True)
            self.lbl_load_val.setStyleSheet("color: black;"); self.lbl_loss_val.setStyleSheet("color: black;")
        else:
            self.btn_mode.setText("ADVANCED MODE: OFF")
            self.btn_mode.setStyleSheet("background-color: #E6E6E6; font-weight: bold; padding: 5px;")
            self.slider_load.setEnabled(False); self.slider_loss.setEnabled(False)
            self.lbl_load_val.setStyleSheet("color: gray;"); self.lbl_loss_val.setStyleSheet("color: gray;")
        self.update_sim()

    def reset_simulation(self):
        self.slider_p.setValue(750)
        self.slider_pf.setValue(int(np.degrees(np.arccos(0.9))))
        self.slider_xd.setValue(128)
        if self.btn_mode.isChecked(): self.btn_mode.click()
        self.slider_load.setValue(100)
        self.slider_loss.setValue(0)
        if hasattr(self.toolbar, '_nav_stack'): self.toolbar._nav_stack.clear()
        self.toolbar.push_current()

    def on_scroll(self, event):
        if event.inaxes != self.ax_phasor: return
        ax = self.ax_phasor
        cur_xlim = ax.get_xlim(); cur_ylim = ax.get_ylim()
        xdata = event.xdata; ydata = event.ydata
        if xdata is None: return
        scale = 1/1.2 if event.button == 'up' else 1.2
        new_w = (cur_xlim[1]-cur_xlim[0])*scale; new_h = (cur_ylim[1]-cur_ylim[0])*scale
        relx = (cur_xlim[1]-xdata)/(cur_xlim[1]-cur_xlim[0])
        rely = (cur_ylim[1]-ydata)/(cur_ylim[1]-cur_ylim[0])
        ax.set_xlim([xdata-new_w*(1-relx), xdata+new_w*relx])
        ax.set_ylim([ydata-new_h*(1-rely), ydata+new_h*rely])
        self.canvas.draw()

    def update_sim(self):
        P_MW = self.slider_p.value()
        PF_Deg = self.slider_pf.value()
        Xd_pu = self.slider_xd.value() / 100.0
        is_adv = self.btn_mode.isChecked()
        load_factor = self.slider_load.value() / 100.0 if is_adv else 1.0
        loss_adder = self.slider_loss.value() / 1000.0 if is_adv else 0.0
        
        self.lbl_p_val.setText(f"{P_MW} MW"); self.lbl_xd_val.setText(f"j{Xd_pu:.2f}")
        self.lbl_load_val.setText(f"{load_factor:.2f}x"); self.lbl_loss_val.setText(f"+{loss_adder:.3f}")
        
        phi_rad = np.radians(PF_Deg); PF = np.cos(phi_rad)
        self.lbl_pf_val.setText(f"{PF:.2f} {'Lead' if PF_Deg>0 else 'Lag'}")
        self.lbl_pf_val.setStyleSheet("color: red;" if PF_Deg>0 else "color: blue;")
            
        P_eff = P_MW * load_factor
        X_gen_sys = Xd_pu * (self.S_base_sys / self.S_gen_rated)
        I_mag = (P_eff * 1e6) / (self.S_base_sys * max(PF, 0.1))
        I_phasor = I_mag * (np.cos(phi_rad) + 1j * np.sin(phi_rad))
        
        V_grid = 1.0 + 0j
        Z_total = complex(self.Z_xfm_pu.real + loss_adder, self.Z_xfm_pu.imag)
        V_term = V_grid + I_phasor * Z_total
        E_int = V_term + I_phasor * 1j * X_gen_sys
        
        # Gerilim Düşümü Bileşenleri (Grafik için)
        V_drop_IR = I_phasor * Z_total.real
        V_drop_IX = I_phasor * 1j * Z_total.imag
        V_drop_Gen = I_phasor * 1j * X_gen_sys

        abs_Eaf = abs(E_int); abs_Vg = abs(V_term)
        vr = ((abs_Eaf - abs_Vg) / abs_Vg) * 100
        p_loss = (I_mag**2) * Z_total.real
        p_out = (P_eff * 1e6) / self.S_base_sys
        eff = (p_out / (p_out + p_loss) * 100) if p_out > 0 else 0
            
        self.lbl_vr.setText(f"Voltage Reg: {vr:.2f} %")
        self.lbl_eff.setText(f"Efficiency: {eff:.2f} % (Loss: {p_loss:.4f} pu)")
        
        self.plot_phasor(V_grid, V_drop_IR, V_drop_IX, V_drop_Gen, E_int, I_phasor, V_term)
        self.plot_bar(V_term, E_int, I_mag)

    def plot_phasor(self, V_grid, V_drop_IR, V_drop_IX, V_drop_Gen, E_int, I_phasor, V_term):
        ax = self.ax_phasor; ax.clear(); ax.grid(True, linestyle='--')
        ax.set_title("Dynamic Phasor Diagram")
        
        def draw_arrow(start, vec, col, lbl=None):
            ax.quiver(start.real, start.imag, vec.real, vec.imag, angles='xy', scale_units='xy', scale=1, color=col, label=lbl, width=0.005, headwidth=4)
        
        origin = 0+0j
        
        # Ana Vektörler
        draw_arrow(origin, V_grid, 'black', 'Grid')
        draw_arrow(origin, V_term, 'blue', 'Vt')
        draw_arrow(origin, E_int, 'cyan', 'Eaf')
        
        # Düşüm Vektörleri
        draw_arrow(V_grid, V_drop_IR, '#D9B310')
        draw_arrow(V_grid + V_drop_IR, V_drop_IX, '#D9541A')
        draw_arrow(V_term, V_drop_Gen, 'magenta')
        draw_arrow(origin, I_phasor*0.4, 'green', 'I')
        
        # Ana Etiketler
        ax.text(V_grid.real/2, -0.05, '$V_{grid}$', ha='center', va='top')
        ax.text(V_term.real/2, V_term.imag/2+0.05, '$V_t$', color='blue')
        ax.text(E_int.real/2, E_int.imag/2+0.2, '$E_{af}$', color='cyan')
        
        # -- YENİ: ARA BİLEŞEN ETİKETLERİ --
        mid_IR = V_grid + V_drop_IR/2
        ax.text(mid_IR.real, mid_IR.imag - 0.10, '$IR$', color='#D9B310', fontsize=9, fontweight='bold')
        
        mid_IX = V_grid + V_drop_IR + V_drop_IX/2
        ax.text(mid_IX.real + 0.05, mid_IX.imag, '$jIX_T$', color='#D9541A', fontsize=9, fontweight='bold')
        
        mid_IXg = V_term + V_drop_Gen/2
        ax.text(mid_IXg.real + 0.05, mid_IXg.imag + 0.02, '$jIX_g$', color='magenta', fontsize=9, fontweight='bold')

        # -- YENİ: CURRENT (I) ETİKETİ EKLENDİ --
        I_end = I_phasor * 0.4
        ax.text(I_end.real, I_end.imag + 0.05, '$I$', color='green', fontsize=10, fontweight='bold')

        m = max(abs(E_int), 1.2) + 0.2
        ax.set_xlim(-0.1, m); ax.set_ylim(-0.2, m)
        ax.legend(loc='upper left', fontsize='small')
        self.canvas.draw()
        if hasattr(self.toolbar, '_nav_stack'): self.toolbar.update()

    def plot_bar(self, V_term, E_int, I_mag):
        ax = self.ax_bar
        ax.clear()
        
        # Değerler
        v_vals = [abs(V_term), 1.0, abs(E_int), I_mag]
        labels = ['$V_t$', '$V_{grid}$', '$E_{af}$', '$I$']
        x_pos = [1, 2, 3, 4]
        
        # Varsayılan Renkler (Normal Durum)
        # Vt: Mavi, Grid: Siyah, Eaf: Kırmızı, I: Yeşil
        bar_colors = ['blue', 'black', 'cyan', 'green']
        title = "System Values (pu)"
        
        # KONTROL: Eğer Under-Excited ise (Eaf < Vt)
        if abs(E_int) < abs(V_term): 
            title = "RISK: Under-Excited!"
            # İsteğin üzerine HEPSİNİ kırmızı yapıyoruz
            bar_colors = ['red', 'red', 'red', 'red']
            
        # Grafiği çiz
        bars = ax.bar(x_pos, v_vals, color=bar_colors)
        
        ax.set_xticks(x_pos)
        ax.set_xticklabels(labels)
        ax.set_title(title)
        ax.set_ylim(0, max(v_vals) * 1.2)
        
        for bar in bars:
            h = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., 1.01*h, f'{h:.3f}', ha='center', va='bottom')
        
        self.canvas.draw()

    # =================================================================
    # --- INFO SEKMESİ VE ŞEMA ÇİZİMİ ---
    # =================================================================
    def setup_info_tab(self):
        layout = QVBoxLayout(self.tab_info)
        
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_content = QWidget()
        scroll_layout = QVBoxLayout(scroll_content)

        title = QLabel("EXPLANATIONS OF THE SYSTEM VARIABLES")
        title.setFont(QFont("Arial", 14, QFont.Weight.Bold))
        title.setStyleSheet("color: #000000; margin-bottom: 10px;")
        title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        scroll_layout.addWidget(title)

        info_texts = [
            "<b>Per-Unit (pu):</b> A method used to analyze power system components with different power (800 MVA and 850 MVA) and voltage (26 kV and 345 kV) ratings on a common mathematical scale. In this study, all system parameters have been normalized using a 850 MVA system base.",
            "<b>Efficiency (η):</b> The ratio of useful power output to the total power generated by the system.A dynamic performance metric integrated into the software's Advanced Mode to evaluate the ratio of useful power output to total generation. Rather than being a fixed value, it is calculated based on real-time system losses-such as the transformer's I<sup>2</sup>R loss.",
            "<b>Voltage Regulation:</b> The percentage difference between the operating voltage under full load and the voltage under no-load conditions. It is a critical parameter for measuring the voltage stability of the generator against load variations.",
            "<b>Transformer Impedance (Z<sub>T</sub>):</b> The total electrical opposition caused by the combined resistance and reactance of the transformer windings. Defined as 0.0025 + j0.057 per-unit in this problem, this value leads to voltage drops and power losses across the system.",
            "<b>Generator Internal Voltage (E<sub>int</sub>):</b> The pure electromotive force induced in the generator windings by the magnetic field. This value represents the voltage behind the generator reactance and was calculated as 34.14 kV (1.313 pu) under the capacitive load condition.",
            "<b>Generator Terminal Voltage (V<sub>t</sub>):</b> The voltage measured at the generator output terminals, just before the transformer input. Following the system analysis, this value was determined to be 25.48 kV (0.980 pu).",
            "<b>Generator Reactance (X<sub>g</sub>):</b> The magnetic opposition to alternating current inherent in the internal structure of the generator. Originally specified as j1.28 per-unit on the generator base, this value was converted to 1.36 per-unit on the system base for the analysis.",
        ]

        for text in info_texts:
            lbl = QLabel(text); lbl.setFont(QFont("Arial", 10)); lbl.setWordWrap(True)
            scroll_layout.addWidget(lbl)
        
        scroll_area.setWidget(scroll_content)
        scroll_area.setMaximumHeight(500)
        layout.addWidget(scroll_area)

        sep = QLabel("SINGLE PHASE CIRCUIT DIAGRAM OF THE SYSTEM (PER-UNIT)")
        sep.setFont(QFont("Arial", 12, QFont.Weight.Bold))
        sep.setStyleSheet("color: #000000; margin-top: 10px; margin-bottom: 0px;")
        sep.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(sep)
        sep.raise_()

        self.schematic_fig = Figure(figsize=(8, 4), dpi=100)
        self.schematic_canvas = FigureCanvas(self.schematic_fig)
        self.schematic_canvas.setMinimumHeight(100)
        
        layout.addWidget(self.schematic_canvas)
        self.schematic_canvas.lower()
        self.draw_schematic()

    def draw_schematic(self):
        ax = self.schematic_fig.add_subplot(111)
        ax.clear(); ax.set_axis_off(); ax.set_aspect('equal')
        lw = 2; color = 'black'

        def draw_inductor(ax, x_start, x_end, y, coils=4, amp=0.2):
            x = np.linspace(x_start, x_end, 100)
            y_coil = y + amp * np.sin(np.linspace(0, coils * 2 * np.pi, 100))
            ax.plot(x, y_coil, color=color, lw=lw)

        def draw_vertical_inductor(ax, x, y_start, y_end, coils=3, amp=0.2):
            y = np.linspace(y_start, y_end, 100)
            x_coil = x + amp * np.sin(np.linspace(0, coils * 2 * np.pi, 100))
            ax.plot(x_coil, y, color=color, lw=lw)

        def draw_resistor(ax, x_start, x_end, y, zigzags=3, amp=0.25):
            x = np.linspace(x_start, x_end, 100)
            y_res = y + amp * np.arcsin(np.sin(np.linspace(0, zigzags * 2 * np.pi, 100))) * (2/np.pi)
            ax.plot(x, y_res, color=color, lw=lw)

        # --- ŞEMA DÜZENİ ---
        ax.plot([1, 9.5], [0, 0], color=color, lw=lw)
        ax.plot([1, 1], [1.4, 2], color=color, lw=lw)
        ax.plot([1, 1], [3, 3.5], color=color, lw=lw)
        ax.plot([1, 3], [3.5, 3.5], color=color, lw=lw)
        ax.plot([7, 9.5], [3.5, 3.5], color=color, lw=lw)
        ax.plot([9.5, 9.5], [3.5, 3.5], color=color, lw=lw) 
        ax.plot([9.5, 9.5], [0.5, 0], color=color, lw=lw)
        ax.plot([1.0, 1.0], [0.6, 0], color=color, lw=lw)
        ax.plot([9.5, 9.5], [3.5, 3], color=color, lw=lw)
        ax.plot([7, 6.8], [3.5, 3.5], color=color, lw=lw)
        ax.plot([3.2, 3], [3.5, 3.5], color=color, lw=lw)

        # 1. Jeneratör
        gen_circle = mpatches.Circle((1, 1), 0.4, fill=False, color=color, lw=lw)
        ax.add_patch(gen_circle)
        t = np.linspace(0.7, 1.3, 50)
        ax.plot(t, 1 + 0.15*np.sin((t-0.7)*4*np.pi), color=color, lw=lw-0.5)
        ax.text(0.4, 1, "$E_{int}$", fontsize=12, fontweight='bold', va='center', ha='right')
        
        # 2. Bobin 
        draw_vertical_inductor(ax, 1, 2, 3, coils=3)
        ax.text(0.4, 2.5, "1.36 pu $X_g$", fontsize=11, ha='right', va='center')
        
        # 3. Vi Oku ve Dot
        ax.plot(1.5, 4.5, marker='o', markersize=6, color='white') 
        ax.annotate("", xy=(1.8, 0), xytext=(1.8, 3.5), arrowprops=dict(arrowstyle="<->", lw=lw, color='black'))
        ax.text(2, 1.75, "$V_t$", fontsize=14, fontweight='bold', va='center')

        # 4. Trafo
        box = mpatches.Rectangle((3, 3), 4, 1, fill=False, color=color, lw=lw)
        ax.add_patch(box)
        draw_resistor(ax, 3.2, 4.8, 3.5)
        draw_inductor(ax, 5.2, 6.8, 3.5)
        ax.plot([4.8, 5.2], [3.5, 3.5], color=color, lw=lw)
        ax.text(5, 4.2, "0.0025+j0.057 pu", fontsize=11, ha='center') 
        ax.text(5, 2.8, "Transformer\nImpedance\n($Z_T$)", fontsize=11, ha='center', va='top')

        # 5. Yük
        load_box = mpatches.Rectangle((9, 0.5), 1, 2.5, fill=False, color=color, lw=lw)
        ax.add_patch(load_box)
        ax.text(9.5, 1.75, "L\nO\nA\nD", fontsize=12, fontweight='bold', ha='center', va='center', linespacing=1.8)

        ax.set_xlim(-0.5, 11); ax.set_ylim(-0.5, 5)
        self.schematic_canvas.draw()

    def setup_about_tab(self):
        self.tab_about_layout = QVBoxLayout(self.tab_about)
        self.anim_container = QWidget()
        self.anim_container.setStyleSheet("background-color: white; border: 1px solid #ddd;")
        self.anim_container.setFixedHeight(750)
        self.tab_about_layout.addWidget(self.anim_container)
        self.scroll_lbl = QLabel(self.anim_container)
        credits_text = (
            "Developed By:\nAli ÇİFT\nEmre YILMAZGÖZ\nBerfin KARABOĞA\n\n"
            "Advisor:\nAli Can ERÜST\n\n"
            "Course:\nElectromechanical Energy Conversion\n\n"
            "Muğla Sıtkı Koçman University\nElectrical & Electronics Engineering\n\n"
            "Project Based On:\nFitzgerald & Kingsley's Electric Machinery\nProblem 2.52 (7th Edition)\n\n"
            "Simulation Details:\nThis software models the interaction between\na synchronous generator and the power grid.\n\n"
            "Features:\n- Dynamic Phasor Diagrams\n- Voltage Regulation Calculation\n- Efficiency Analysis\n\n"
            "Date:\nDecember 2025\n\n© 2025 MSKU Engineering Faculty"
        )
        self.scroll_lbl.setText(credits_text)
        self.scroll_lbl.setAlignment(Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignTop)
        self.scroll_lbl.setFont(QFont("Arial", 14))
        self.scroll_lbl.setStyleSheet("color: #555;")
        self.scroll_lbl.resize(1600, 1200)
        self.timer = QTimer()
        self.timer.timeout.connect(self.animate_text)
        self.timer.start(40)
        self.scroll_y = 400

    def animate_text(self):
        self.scroll_y -= 1.5
        if self.scroll_y < -700: self.scroll_y = self.anim_container.height()
        center_x = (self.anim_container.width() - self.scroll_lbl.width()) // 2
        self.scroll_lbl.move(int(center_x), int(self.scroll_y))

if __name__ == '__main__':
    app = QApplication(sys.argv)
    font = QFont("Segoe UI", 9)
    app.setFont(font)
    window = GeneratorAnalyzer()
    window.show()
    sys.exit(app.exec())