"""
A utility module that will hold the standardized SymPortal color schemes and tools
"""
import random

# These are the top 40 most abundnant named sequences in the SymPortal db.
# They have a color hardcoded to them.
pre_def_color_dict = {
    'A1': "#FFFF00", 'C3': "#1CE6FF", 'C15': "#FF34FF", 'A1bo': "#FF4A46", 'D1': "#008941",
    'C1': "#006FA6", 'C27': "#A30059", 'D4': "#FFDBE5", 'C3u': "#7A4900", 'C42.2': "#0000A6",
    'A1bp': "#63FFAC", 'C115': "#B79762", 'C1b': "#004D43", 'C1d': "#8FB0FF", 'A1c': "#997D87",
    'C66': "#5A0007", 'A1j': "#809693", 'B1': "#FEFFE6", 'A1k': "#1B4400", 'A4': "#4FC601",
    'A1h': "#3B5DFF", 'C50a': "#4A3B53", 'C39': "#FF2F80", 'C3dc': "#61615A", 'D4c': "#BA0900",
    'C3z': "#6B7900", 'C21': "#00C2A0", 'C116': "#FFAA92", 'A1cc': "#FF90C9", 'C72': "#B903AA",
    'C15cl': "#D16100", 'C31': "#DDEFFF", 'C15cw': "#000035", 'A1bv': "#7B4F4B", 'D6': "#A1C299",
    'A4m': "#300018", 'C42a': "#0AA6D8", 'C15cr': "#013349", 'C50l': "#00846F", 'C42g': "#372101"}


# These are the colors used by SymPortal to represent sequences that are not in the pre_def_color_dict
# This list does not contain colors used in the pre_def_color_dict
color_list = [
    "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09", "#00489C", "#6F0062",
    "#0CBD66",
    "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66", "#885578", "#FAD09F", "#FF8A9A", "#D157A0",
    "#BEC459",
    "#456648", "#0086ED", "#886F4C", "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9",
    "#FF913F",
    "#938A81", "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
    "#7900D7",
    "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700", "#549E79", "#FFF69F",
    "#201625",
    "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329", "#5B4534", "#FDE8DC", "#404E55", "#0089A3",
    "#CB7E98",
    "#A4E804", "#324E72", "#6A3A4C", "#83AB58", "#001C1E", "#D1F7CE", "#004B28", "#C8D0F6", "#A3A489",
    "#806C66",
    "#222800", "#BF5650", "#E83000", "#66796D", "#DA007C", "#FF1A59", "#8ADBB4", "#1E0200", "#5B4E51",
    "#C895C5",
    "#320033", "#FF6832", "#66E1D3", "#CFCDAC", "#D0AC94", "#7ED379", "#012C58", "#7A7BFF", "#D68E01",
    "#353339",
    "#78AFA1", "#FEB2C6", "#75797C", "#837393", "#943A4D", "#B5F4FF", "#D2DCD5", "#9556BD", "#6A714A",
    "#001325",
    "#02525F", "#0AA3F7", "#E98176", "#DBD5DD", "#5EBCD1", "#3D4F44", "#7E6405", "#02684E", "#962B75",
    "#8D8546",
    "#9695C5", "#E773CE", "#D86A78", "#3E89BE", "#CA834E", "#518A87", "#5B113C", "#55813B", "#E704C4",
    "#00005F",
    "#A97399", "#4B8160", "#59738A", "#FF5DA7", "#F7C9BF", "#643127", "#513A01", "#6B94AA", "#51A058",
    "#A45B02",
    "#1D1702", "#E20027", "#E7AB63", "#4C6001", "#9C6966", "#64547B", "#97979E", "#006A66", "#391406",
    "#F4D749",
    "#0045D2", "#006C31", "#DDB6D0", "#7C6571", "#9FB2A4", "#00D891", "#15A08A", "#BC65E9", "#FFFFFE",
    "#C6DC99",
    "#203B3C", "#671190", "#6B3A64", "#F5E1FF", "#FFA0F2", "#CCAA35", "#374527", "#8BB400", "#797868",
    "#C6005A",
    "#3B000A", "#C86240", "#29607C", "#402334", "#7D5A44", "#CCB87C", "#B88183", "#AA5199", "#B5D6C3",
    "#A38469",
    "#9F94F0", "#A74571", "#B894A6", "#71BB8C", "#00B433", "#789EC9", "#6D80BA", "#953F00", "#5EFF03",
    "#E4FFFC",
    "#1BE177", "#BCB1E5", "#76912F", "#003109", "#0060CD", "#D20096", "#895563", "#29201D", "#5B3213",
    "#A76F42",
    "#89412E", "#1A3A2A", "#494B5A", "#A88C85", "#F4ABAA", "#A3F3AB", "#00C6C8", "#EA8B66", "#958A9F",
    "#BDC9D2",
    "#9FA064", "#BE4700", "#658188", "#83A485", "#453C23", "#47675D", "#3A3F00", "#061203", "#DFFB71",
    "#868E7E",
    "#98D058", "#6C8F7D", "#D7BFC2", "#3C3E6E", "#D83D66", "#2F5D9B", "#6C5E46", "#D25B88", "#5B656C",
    "#00B57F",
    "#545C46", "#866097", "#365D25", "#252F99", "#00CCFF", "#674E60", "#FC009C", "#92896B"]

# These are the standard greys used in SymPortal figure plotting for sequences
greys = ['#D0CFD4', '#89888D', '#4A4A4C', '#8A8C82', '#D4D5D0', '#53544F']

genus_color_dict = {
                    # Pastel red
                    'A': "#FF7F7F",
                    # Pastel yellow
                    'B': "#FFFF7F",
                    # Pastel green
                    'C': "#BFFF7F",
                    # Pastel Blue
                    'D': "#7F7FFF",
                    # Pastel Purple
                    'E': "#BF7FFF",
                    # Pastel orange
                    'F': "#FFBF7F",
                    # Pastel Cyan
                    'G': "#7FFFFF",
                    # Pastel Magenta
                    'H': "#FF7FFF",
                    # Grey
                    'I': "#74737A"
                }

def create_color_list(
        sq_dist_cutoff=None, mix_col=None, num_cols=50, time_out_iterations=10000,
        avoid_black_and_white=True, error_on_fail=False, warnings_off=False
):
    """
    Returns a random list of colors in a format that can be used for plotting in matplotlib.

    :param sq_dist_cutoff: A minimum distance between colors. [None]

    :param mix_col: A color that all generated colors will be 'mixed' with. This is often set to white to make
    pastel colors for representing ITS2 type profile colors. [None]

    :param num_cols: The number of colors to be generated [50]

    :param time_out_iterations: The number of iterations to attempt before either
    raising a RuntimeError (if error_on_fail is True) or reducing the sq_dist_cutoff param by 10%
    and trying to generate the colors again (if error_on_fail is False). [1000]

    :param avoid_black_and_white: Whether to allow white and black in the colors generated.
    If True, all generated colors must be sq_dist_cutoff away from white and black in
    addition to the other colors generated. [True]

    :param error_on_fail: Whether to raise a RuntimeError if the color list cannot be generated using the above params.
    see time_out_iterations. [False]

    :param warnings_off: Whether to raise RuntimeWarning when reducing the number of
    time_out_iterations after failing to generate the requested number of colors.
    If True, no RuntimeWarning will be raised. [False]

    :return: A list of colors as tuples of R, G, B values that can be used for plotting in matplotlib.

    Example usage:
    color_palette_pas_gen = (
                                '#%02x%02x%02x' % rgb_tup for rgb_tup in
                                self._create_color_list(
                                    mix_col=(255, 255, 255),
                                    sq_dist_cutoff=5000,
                                    num_cols=8,
                                    time_out_iterations=10000
                                )
                            )
    """
    new_colors = []
    min_dist = []
    attempt = 0
    while len(new_colors) < num_cols:
        attempt += 1
        # Check to see if we have run out of iteration attempts to find a color that fits into the color space
        if attempt > time_out_iterations:
            if error_on_fail:
                raise RuntimeError(
                    f'Colour generation timed out. We have tried {attempt} iterations of color generation '
                    'and have not been able to find a color that fits into your defined color space.\n'
                    'Please lower the number of colors you are trying to find, '
                    'the minimum distance between them, or both.')
            else:
                # Reduce the sq_dist_cutoff by 10% and reset attempts and continue to try
                old_sq_dist_cutoff = sq_dist_cutoff
                sq_dist_cutoff = int(sq_dist_cutoff - (sq_dist_cutoff * 0.1))
                attempt = 0
                if not warnings_off:
                    raise RuntimeWarning(
                        f'We have tried {attempt} iterations of color generation '
                        f'and have not been able to find a color that fits into your defined color space.\n'
                        f'Automatically lowering your distance from {old_sq_dist_cutoff} to {sq_dist_cutoff}.\n'
                        f'Continuing with color generation'
                    )
        if mix_col:
            r = int((random.randint(0, 255) + mix_col[0]) / 2)
            g = int((random.randint(0, 255) + mix_col[1]) / 2)
            b = int((random.randint(0, 255) + mix_col[2]) / 2)
        else:
            r = random.randint(0, 255)
            g = random.randint(0, 255)
            b = random.randint(0, 255)

        # now check to see whether the new color is within a given distance
        # if the avoids are true also
        good_dist = True
        if sq_dist_cutoff:
            dist_list = []
            for i in range(len(new_colors)):
                distance = (new_colors[i][0] - r) ** 2 + (new_colors[i][1] - g) ** 2 + (
                        new_colors[i][2] - b) ** 2
                dist_list.append(distance)
                if distance < sq_dist_cutoff:
                    good_dist = False
                    break
            # now check against black and white
            d_to_black = (r - 0) ** 2 + (g - 0) ** 2 + (b - 0) ** 2
            d_to_white = (r - 255) ** 2 + (g - 255) ** 2 + (b - 255) ** 2
            if avoid_black_and_white:
                if d_to_black < sq_dist_cutoff or d_to_white < sq_dist_cutoff:
                    good_dist = False
            if dist_list:
                min_dist.append(min(dist_list))
        if good_dist:
            new_colors.append((r, g, b))
            attempt = 0

    return new_colors
