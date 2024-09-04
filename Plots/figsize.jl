## helperfunctions for setting correct size and fontsize of plot

function size_pt(; width = 455.24411)

    if width == :thesis
        width_pt = 426.79135
    elseif width == :beamer
        width_pt = 307.28987
    else
        width_pt = width
    end

    inches_per_pt = 1 / 72.27

    golden_ratio = (sqrt(5) - 1) / 2

    width_in = width_pt * inches_per_pt

    height_in = width_in * golden_ratio

    return (width_in, height_in)
end