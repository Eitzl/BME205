"""
Power iteration on A = [[-5, 2], [-9, 6]], animated with Manim.

Render:  manim -qh power_iteration.py PowerIteration
"""
import numpy as np
from manim import *

A = np.array([[-5.0, 2.0], [-9.0, 6.0]])
E1 = np.array([2.0, 9.0]) / np.sqrt(85)   # eigenvector, lambda = 4  (dominant)
E2 = np.array([1.0, 1.0]) / np.sqrt(2)    # eigenvector, lambda = -3

GOLD_C = "#F5B83D"
BLUE_C = "#4FA3F7"
BG = "#0E1117"


def to3(v):
    return np.array([v[0], v[1], 0.0])


class PowerIteration(Scene):
    def construct(self):
        self.camera.background_color = BG
        self.intro()
        self.swarm()
        self.explain()
        self.outro()

    # ------------------------------------------------------------------ 1
    def intro(self):
        title = Text("Power Iteration", weight=BOLD).scale(1.1)
        sub = Text("why a random vector finds the eigenvector", color=GREY_B).scale(0.5)
        sub.next_to(title, DOWN)
        self.play(Write(title), FadeIn(sub, shift=UP * 0.2))
        self.wait(0.5)
        self.play(VGroup(title, sub).animate.to_edge(UP))

        mat = MathTex(r"A = \begin{bmatrix} -5 & 2 \\ -9 & 6 \end{bmatrix}").scale(1.1)
        recipe = VGroup(
            MathTex(r"1.\ \ v_0 = \text{random vector}"),
            MathTex(r"2.\ \ v_{k+1} = \frac{A v_k}{\lVert A v_k \rVert}"),
            MathTex(r"3.\ \ \text{repeat}"),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.35).scale(0.85)
        VGroup(mat, recipe).arrange(RIGHT, buff=1.4).shift(DOWN * 0.3)
        self.play(FadeIn(mat))
        self.play(LaggedStart(*[Write(r) for r in recipe], lag_ratio=0.5))
        q = Text("Where does it end up?", color=YELLOW).scale(0.6).to_edge(DOWN, buff=0.8)
        self.play(FadeIn(q))
        self.wait(1.5)
        self.play(*[FadeOut(m) for m in self.mobjects])

    # ------------------------------------------------------------------ 2
    def swarm(self):
        plane = NumberPlane(
            x_range=[-8, 8, 1], y_range=[-5, 5, 1],
            background_line_style={"stroke_color": BLUE_E, "stroke_opacity": 0.35, "stroke_width": 1},
            axis_config={"stroke_opacity": 0.6},
        )
        R = 2.8
        circle = Circle(radius=R, color=GREY_C, stroke_width=1.5, stroke_opacity=0.6)
        self.play(Create(plane, run_time=1.2), Create(circle))

        # 24 random-looking starting directions, avoiding the exact (1,1) line
        rng = np.random.default_rng(7)
        angles = rng.uniform(0, 2 * PI, 200)
        # keep starts at least ~14 degrees off the (1,1) line so all converge on screen
        off = np.abs(((angles - PI / 4) + PI / 2) % PI - PI / 2)
        angles = np.sort(angles[off > np.radians(14)][:24])
        trackers = [ValueTracker(a) for a in angles]
        colors = color_gradient([PURPLE_A, TEAL_C, PINK, BLUE_B], len(angles))

        def make_arrow(t, c):
            return always_redraw(lambda: Arrow(
                ORIGIN, R * np.array([np.cos(t.get_value()), np.sin(t.get_value()), 0]),
                buff=0, color=c, stroke_width=3, max_tip_length_to_length_ratio=0.08,
            ))

        arrows = VGroup(*[make_arrow(t, c) for t, c in zip(trackers, colors)])
        tips = VGroup(*[
            always_redraw(lambda t=t, c=c: Dot(
                R * np.array([np.cos(t.get_value()), np.sin(t.get_value()), 0]),
                radius=0.05, color=c))
            for t, c in zip(trackers, colors)
        ])
        trails = VGroup(*[
            TracedPath(d.get_center, stroke_color=c, stroke_width=2.5, dissipating_time=0.6)
            for d, c in zip(tips, colors)
        ])
        label = Text("24 random starting vectors", color=GREY_A).scale(0.5)
        label.to_corner(UL).add_background_rectangle(color=BG, opacity=0.8)
        self.add(trails)
        self.play(LaggedStart(*[GrowArrow(a) for a in arrows], lag_ratio=0.04), FadeIn(tips), FadeIn(label))
        self.wait(0.5)

        counter_k = ValueTracker(0)
        counter = always_redraw(lambda: MathTex(
            rf"k = {int(counter_k.get_value())}").scale(0.9).to_corner(UR)
            .add_background_rectangle(color=BG, opacity=0.8))
        self.play(FadeIn(counter), FadeOut(label))

        n_steps = 16
        for k in range(n_steps):
            anims = []
            for t in trackers:
                th = t.get_value()
                v = A @ np.array([np.cos(th), np.sin(th)])
                new = np.arctan2(v[1], v[0])
                d = (new - th + PI) % (2 * PI) - PI    # shortest way round the circle
                anims.append(t.animate.set_value(th + d))
            rt = max(1.1 - 0.07 * k, 0.4)
            self.play(*anims, counter_k.animate.set_value(k + 1), run_time=rt, rate_func=smooth)
            self.wait(0.1)

        # reveal the line they collapsed onto
        e1_line = DashedLine(to3(-4.6 * E1), to3(4.6 * E1), color=GOLD_C, stroke_width=4)
        e1_lab = MathTex(r"\text{direction } (2,\,9)", color=GOLD_C).scale(0.7)
        e1_lab.next_to(to3(3.3 * E1), RIGHT, buff=0.25).add_background_rectangle(color=BG, opacity=0.85)
        self.play(Create(e1_line), FadeIn(e1_lab))
        verdict = MathTex(r"v_k \to \pm \tfrac{1}{\sqrt{85}}(2,\,9) \approx \pm(0.217,\ 0.976)")
        verdict.scale(0.8).to_edge(DOWN, buff=0.4).add_background_rectangle(color=BG, opacity=0.9)
        self.play(Write(verdict))
        self.wait(1.2)
        why = Text("Every start ends up here. Why?", color=YELLOW).scale(0.6)
        why.next_to(verdict, UP, buff=0.25).add_background_rectangle(color=BG, opacity=0.9)
        self.play(FadeIn(why))
        self.wait(1.5)
        self.remove(trails)
        self.play(*[FadeOut(m) for m in self.mobjects if m is not plane], run_time=0.8)
        self.plane = plane

    # ------------------------------------------------------------------ 3
    def explain(self):
        plane = self.plane

        # The two special directions
        e1_line = DashedLine(to3(-4.6 * E1), to3(4.6 * E1), color=GOLD_C, stroke_opacity=0.6)
        e2_line = DashedLine(to3(-6 * E2), to3(6 * E2), color=BLUE_C, stroke_opacity=0.6)
        self.play(Create(e1_line), Create(e2_line))

        facts = VGroup(
            MathTex(r"A\begin{bmatrix}2\\9\end{bmatrix} = ", r"4", r"\begin{bmatrix}2\\9\end{bmatrix}"),
            MathTex(r"A\begin{bmatrix}1\\1\end{bmatrix} = ", r"-3", r"\begin{bmatrix}1\\1\end{bmatrix}"),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.4).scale(0.7)
        facts[0].set_color(GOLD_C)
        facts[1].set_color(BLUE_C)
        facts.to_corner(UL).shift(DOWN * 0.1)
        bg = BackgroundRectangle(facts, color=BG, fill_opacity=0.9, buff=0.2)
        head = Text("Two eigenvectors: A only stretches them", color=GREY_A).scale(0.45)
        head.next_to(bg, UP, buff=0.1, aligned_edge=LEFT)
        hb = BackgroundRectangle(head, color=BG, fill_opacity=0.9, buff=0.08)
        self.play(FadeIn(bg), FadeIn(hb), FadeIn(head), Write(facts))
        self.wait(2)

        # Any vector = gold part + blue part
        a = 1.5
        b = ValueTracker(2.6)
        gold = Arrow(ORIGIN, to3(a * E1), buff=0, color=GOLD_C, stroke_width=6,
                     max_tip_length_to_length_ratio=0.15)
        blue = always_redraw(lambda: Arrow(
            to3(a * E1), to3(a * E1 + b.get_value() * E2), buff=0, color=BLUE_C,
            stroke_width=6, max_tip_length_to_length_ratio=0.12 if abs(b.get_value()) > 0.4 else 0.35))
        v = always_redraw(lambda: Arrow(
            ORIGIN, to3(a * E1 + b.get_value() * E2), buff=0, color=WHITE, stroke_width=5,
            max_tip_length_to_length_ratio=0.1))
        v_lab = always_redraw(lambda: MathTex("v", color=WHITE).scale(0.8)
                              .next_to(v.get_end(), RIGHT if b.get_value() > 0 else LEFT, buff=0.12))

        cap = Text("Split any vector into a gold part + a blue part", color=GREY_A).scale(0.5)
        cap.to_edge(DOWN, buff=0.45).add_background_rectangle(color=BG, opacity=0.9)
        self.play(GrowArrow(v), FadeIn(v_lab), FadeIn(cap))
        self.wait(0.6)
        self.play(GrowArrow(gold), Create(blue))
        self.wait(1.5)

        # The trick: divide by 4 instead of normalizing (same direction)
        cap2 = MathTex(
            r"\text{Each step: } \tfrac{1}{4}A \text{ keeps gold }(\times 1)"
            r"\text{, blue gets } \times(-\tfrac{3}{4})").scale(0.75)
        cap2.to_edge(DOWN, buff=0.45).add_background_rectangle(color=BG, opacity=0.9)
        note = Text("(dividing by 4 instead of normalizing — same direction, easier to see)",
                    color=GREY_B).scale(0.33)
        note.next_to(cap2, UP, buff=0.12).add_background_rectangle(color=BG, opacity=0.9)
        self.play(FadeOut(cap), FadeIn(cap2), FadeIn(note))

        k_tr = ValueTracker(0)
        ratio = always_redraw(lambda: MathTex(
            r"\frac{|\text{blue}|}{|\text{gold}|} \propto \left(\tfrac{3}{4}\right)^{"
            + str(int(k_tr.get_value())) + r"} = "
            + f"{0.75 ** int(k_tr.get_value()):.3f}").scale(0.75).to_corner(UR)
            .add_background_rectangle(color=BG, opacity=0.9))
        self.play(FadeIn(ratio))

        trail = TracedPath(lambda: v.get_end(), stroke_color=WHITE, stroke_width=2,
                           stroke_opacity=0.5)
        self.add(trail)
        for k in range(12):
            rt = max(1.3 - 0.09 * k, 0.45)
            self.play(b.animate.set_value(b.get_value() * -0.75),
                      k_tr.animate.set_value(k + 1), run_time=rt)
            self.wait(0.15)
        self.wait(1)

        punch = Text("The blue part shrinks away. Only the gold direction survives.",
                     color=YELLOW).scale(0.48)
        punch.to_edge(DOWN, buff=0.45).add_background_rectangle(color=BG, opacity=0.95)
        self.play(FadeOut(cap2), FadeOut(note), FadeIn(punch))
        self.wait(2.5)
        self.remove(trail)
        self.play(*[FadeOut(m) for m in self.mobjects])

    # ------------------------------------------------------------------ 4
    def outro(self):
        lines = VGroup(
            Text("Takeaways", weight=BOLD).scale(0.8),
            MathTex(r"\bullet\ \ \text{Power iteration finds the eigenvector with the largest } |\lambda|"),
            MathTex(r"\bullet\ \ \text{Error shrinks like } \left|\lambda_2/\lambda_1\right|^k = (3/4)^k"),
            MathTex(r"\bullet\ \ \lambda_2 = -3 < 0 \text{ makes it zig-zag across the line}"),
            MathTex(r"\bullet\ \ \text{Normalizing only stops the numbers exploding}"),
        )
        for m in lines[1:]:
            m.scale(0.7)
        lines.arrange(DOWN, aligned_edge=LEFT, buff=0.4).move_to(ORIGIN)
        self.play(LaggedStart(*[FadeIn(l, shift=RIGHT * 0.2) for l in lines], lag_ratio=0.35))
        self.wait(3)
