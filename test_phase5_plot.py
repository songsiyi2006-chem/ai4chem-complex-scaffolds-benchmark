"""Render the near-zero-yield case without data-coordinate label explosion."""
import tempfile
import unittest
from pathlib import Path
from PIL import Image
from test_phase3_5_rerun import functions


class PlotTests(unittest.TestCase):
    def test_near_zero_yield_bounded_canvas(self):
        names = ['R', 'Cat', 'RC', 'I1Cat', 'I1', 'P_R', 'P_S', 'P_elim', 'Q', 'P_poly']
        with tempfile.TemporaryDirectory() as temp:
            data = {'profile_298K': {'t': [1e-9, 1., 1e5],
                    'y': {name: [0.01]*3 for name in names}},
                    'T_sweep': [{'T': t, 'ee': 0., 'yield': (t-249)*1e-22}
                                for t in range(250, 351, 10)], 'ee_curtin_298K': 0.}
            ns = functions(5, ['figure_2'], RESULTS={}, XB=data, FIG=Path(temp),
                           style_axes=lambda ax: None, _log=lambda msg: None)
            ns['figure_2']()
            with Image.open(Path(temp)/'fig2_stiff_microkinetics_profile.png') as img:
                self.assertLess(img.width, 6000)
                self.assertLess(img.height, 4000)


if __name__ == '__main__':
    unittest.main()
