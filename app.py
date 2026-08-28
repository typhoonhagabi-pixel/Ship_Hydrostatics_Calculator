from flask import Flask, render_template, request, jsonify

app = Flask(__name__)


# ============================================================
# UNIT CONVERSION
# Internal calculation system:
# SI -> metres, m2, m3, tonnes
# ============================================================

def to_si_length(value, unit_system):
    """Convert length from selected system to metres."""
    if unit_system == "FPS":
        return value * 0.3048
    return value


def from_si_length(value, unit_system):
    """Convert metres to selected length unit."""
    if unit_system == "FPS":
        return value / 0.3048
    return value


def from_si_area(value, unit_system):
    """Convert m2 to selected area unit."""
    if unit_system == "FPS":
        return value / 0.09290304
    return value


def from_si_volume(value, unit_system):
    """Convert m3 to selected volume unit."""
    if unit_system == "FPS":
        return value / 0.028316846592
    return value


def from_si_displacement(value, unit_system):
    """
    SI displacement is tonnes.
    FPS output is long tons.
    """
    if unit_system == "FPS":
        return value / 1.0160469088
    return value


def to_si_density(value, unit_system):
    """
    SI density:
        tonnes / m3

    FPS density:
        lb / ft3
    """
    if unit_system == "FPS":
        # 1 lb/ft3 = 0.016018463 t/m3
        return value * 0.016018463
    return value


def from_si_density(value, unit_system):
    """Convert t/m3 to selected density."""
    if unit_system == "FPS":
        return value / 0.016018463
    return value


def from_si_tpc(value, unit_system):
    """
    SI:
        tonnes/cm

    FPS:
        long tons/in
    """
    if unit_system == "FPS":
        # 1 tonne/cm -> approximately 0.0445 long ton/in
        return value * 0.044505
    return value


# ============================================================
# SIMPSON 1/3 MULTIPLIERS
# ============================================================

def simpson_multipliers(n):

    if n < 3:
        raise ValueError(
            "At least 3 stations are required."
        )

    if n % 2 == 0:
        raise ValueError(
            "Number of stations must be odd."
        )

    multipliers = []

    for i in range(n):

        if i == 0 or i == n - 1:
            multipliers.append(1)

        elif i % 2 == 1:
            multipliers.append(4)

        else:
            multipliers.append(2)

    return multipliers


# ============================================================
# HYDROSTATIC CALCULATIONS
# ============================================================

def calculate_hydrostatics(
    lpp,
    drafts,
    offsets,
    density
):

    number_of_waterlines = len(drafts)
    number_of_stations = len(offsets[0])

    # Distance between stations
    h = lpp / (number_of_stations - 1)

    multipliers = simpson_multipliers(
        number_of_stations
    )

    # --------------------------------------------------------
    # Check waterline drafts
    # --------------------------------------------------------

    for i in range(len(drafts) - 1):

        if drafts[i + 1] <= drafts[i]:

            raise ValueError(
                "Draft values must increase from "
                "one waterline to the next."
            )

    # --------------------------------------------------------
    # WATERPLANE AREA
    # --------------------------------------------------------

    waterplane_areas = []

    moments_of_inertia = []

    lcf_values = []

    for w in range(number_of_waterlines):

        half_breadths = offsets[w]

        # Area
        sum_aw = sum(
            y * m
            for y, m in zip(
                half_breadths,
                multipliers
            )
        )

        aw = (
            2
            * h
            / 3
            * sum_aw
        )

        waterplane_areas.append(aw)

        # ----------------------------------------------------
        # Moment of inertia about centreline
        #
        # I = integral(2/3 y^3 dx)
        # ----------------------------------------------------

        sum_i = sum(
            (y ** 3) * m
            for y, m in zip(
                half_breadths,
                multipliers
            )
        )

        I = (
            2
            * h
            / 3
            * sum_i
        )

        moments_of_inertia.append(I)

        # ----------------------------------------------------
        # LCF
        # ----------------------------------------------------

        weighted_area = 0.0
        weighted_x = 0.0

        for s in range(number_of_stations):

            x = s * h

            width = 2 * half_breadths[s]

            m = multipliers[s]

            weighted_area += width * m

            weighted_x += width * x * m

        if weighted_area > 0:

            lcf = weighted_x / weighted_area

        else:

            lcf = 0.0

        lcf_values.append(lcf)

    # --------------------------------------------------------
    # VOLUME
    #
    # Waterlines may be odd OR even.
    #
    # We integrate the waterplane area curve using
    # trapezoidal integration because the waterline spacing
    # may be arbitrary.
    # --------------------------------------------------------

    volumes = [0.0]

    cumulative_volume = 0.0

    for w in range(1, number_of_waterlines):

        dz = drafts[w] - drafts[w - 1]

        average_area = (
            waterplane_areas[w]
            + waterplane_areas[w - 1]
        ) / 2.0

        slice_volume = (
            average_area * dz
        )

        cumulative_volume += slice_volume

        volumes.append(
            cumulative_volume
        )

    # --------------------------------------------------------
    # HYDROSTATIC TABLE
    # --------------------------------------------------------

    results = []

    for w in range(number_of_waterlines):

        draft = drafts[w]

        volume = volumes[w]

        aw = waterplane_areas[w]

        I = moments_of_inertia[w]

        lcf = lcf_values[w]

        # ----------------------------------------------------
        # Displacement
        # ----------------------------------------------------

        displacement = (
            volume * density
        )

        # ----------------------------------------------------
        # TPC
        #
        # tonnes per centimetre immersion
        # ----------------------------------------------------

        tpc = (
            aw
            * density
            / 100.0
        )

        # ----------------------------------------------------
        # BM
        # ----------------------------------------------------

        if volume > 0:

            bm = I / volume

        else:

            bm = 0.0

        # ----------------------------------------------------
        # KB
        #
        # Vertical centre of buoyancy calculated from
        # volume slices.
        # ----------------------------------------------------

        if volume > 0:

            vertical_moment = 0.0

            for j in range(1, w + 1):

                z1 = drafts[j - 1]
                z2 = drafts[j]

                a1 = waterplane_areas[j - 1]
                a2 = waterplane_areas[j]

                dz = z2 - z1

                slice_volume = (
                    (a1 + a2)
                    / 2.0
                    * dz
                )

                # Centroid approximation
                slice_z = (
                    z1 + z2
                ) / 2.0

                vertical_moment += (
                    slice_volume
                    * slice_z
                )

            kb = (
                vertical_moment
                / volume
            )

        else:

            kb = 0.0

        # ----------------------------------------------------
        # KM
        # ----------------------------------------------------

        km = kb + bm

        # ----------------------------------------------------
        # LCB
        #
        # Calculate longitudinal centre of buoyancy.
        # ----------------------------------------------------

        if volume > 0:

            longitudinal_moment = 0.0

            for s in range(number_of_stations):

                x = s * h

                station_areas = []

                for j in range(w + 1):

                    section_area = (
                        2
                        * offsets[j][s]
                    )

                    station_areas.append(
                        section_area
                    )

                station_volume = 0.0

                for j in range(1, w + 1):

                    dz = (
                        drafts[j]
                        - drafts[j - 1]
                    )

                    station_volume += (
                        (
                            station_areas[j]
                            + station_areas[j - 1]
                        )
                        / 2.0
                        * dz
                    )

                longitudinal_moment += (
                    station_volume * x
                )

            lcb = (
                longitudinal_moment
                / volume
            )

        else:

            lcb = 0.0

        # ----------------------------------------------------
        # Block coefficient
        #
        # Cb = volume / (L x B x T)
        #
        # B is maximum moulded breadth.
        # ----------------------------------------------------

        max_breadth = max(
            offsets[j][s]
            for j in range(w + 1)
            for s in range(number_of_stations)
        ) * 2

        if (
            volume > 0
            and lpp > 0
            and max_breadth > 0
            and draft > 0
        ):

            cb = (
                volume
                /
                (
                    lpp
                    * max_breadth
                    * draft
                )
            )

        else:

            cb = 0.0

        results.append({

            "draft": draft,

            "volume": volume,

            "displacement": displacement,

            "area": aw,

            "tpc": tpc,

            "kb": kb,

            "bm": bm,

            "km": km,

            "lcb": lcb,

            "lcf": lcf,

            "cb": cb

        })

    return results


# ============================================================
# MAIN PAGE
# ============================================================

@app.route("/")
def index():

    return render_template(
        "index.html"
    )


# ============================================================
# CALCULATE API
# ============================================================

@app.route(
    "/calculate",
    methods=["POST"]
)
def calculate():

    try:

        data = request.get_json()

        unit_system = data.get(
            "unit_system",
            "SI"
        )

        if unit_system not in [
            "SI",
            "FPS"
        ]:

            raise ValueError(
                "Invalid unit system."
            )

        vessel_name = data.get(
            "vessel_name",
            "Unnamed Vessel"
        )

        # ----------------------------------------------------
        # Input values
        # ----------------------------------------------------

        lpp_input = float(
            data["lpp"]
        )

        density_input = float(
            data["density"]
        )

        drafts_input = [
            float(x)
            for x in data["drafts"]
        ]

        offsets_input = [

            [
                float(x)
                for x in row
            ]

            for row in data["offsets"]

        ]

        # ----------------------------------------------------
        # Convert inputs to SI
        # ----------------------------------------------------

        lpp = to_si_length(
            lpp_input,
            unit_system
        )

        density = to_si_density(
            density_input,
            unit_system
        )

        drafts = [

            to_si_length(
                x,
                unit_system
            )

            for x in drafts_input

        ]

        offsets = [

            [
                to_si_length(
                    x,
                    unit_system
                )

                for x in row
            ]

            for row in offsets_input

        ]

        # ----------------------------------------------------
        # Basic validation
        # ----------------------------------------------------

        if lpp <= 0:

            raise ValueError(
                "LPP must be greater than zero."
            )

        if density <= 0:

            raise ValueError(
                "Density must be greater than zero."
            )

        if len(drafts) < 2:

            raise ValueError(
                "At least 2 waterlines are required."
            )

        if len(offsets) != len(drafts):

            raise ValueError(
                "Number of waterlines and offset rows "
                "must be the same."
            )

        number_of_stations = len(
            offsets[0]
        )

        # Stations MUST be odd

        if number_of_stations % 2 == 0:

            raise ValueError(
                "Number of stations must be odd."
            )

        # Every waterline must have same station count

        for row in offsets:

            if len(row) != number_of_stations:

                raise ValueError(
                    "Every waterline must contain "
                    "the same number of stations."
                )

        # ----------------------------------------------------
        # Calculate
        # ----------------------------------------------------

        results_si = calculate_hydrostatics(

            lpp=lpp,

            drafts=drafts,

            offsets=offsets,

            density=density

        )

        # ----------------------------------------------------
        # Convert results back to selected system
        # ----------------------------------------------------

        results = []

        for row in results_si:

            results.append({

                "draft":
                    from_si_length(
                        row["draft"],
                        unit_system
                    ),

                "volume":
                    from_si_volume(
                        row["volume"],
                        unit_system
                    ),

                "displacement":
                    from_si_displacement(
                        row["displacement"],
                        unit_system
                    ),

                "area":
                    from_si_area(
                        row["area"],
                        unit_system
                    ),

                "tpc":
                    from_si_tpc(
                        row["tpc"],
                        unit_system
                    ),

                "kb":
                    from_si_length(
                        row["kb"],
                        unit_system
                    ),

                "bm":
                    from_si_length(
                        row["bm"],
                        unit_system
                    ),

                "km":
                    from_si_length(
                        row["km"],
                        unit_system
                    ),

                "lcb":
                    from_si_length(
                        row["lcb"],
                        unit_system
                    ),

                "lcf":
                    from_si_length(
                        row["lcf"],
                        unit_system
                    ),

                "cb":
                    row["cb"]

            })

        # ----------------------------------------------------
        # Station positions
        # ----------------------------------------------------

        x_positions_si = [

            i
            * lpp
            / (number_of_stations - 1)

            for i in range(
                number_of_stations
            )

        ]

        x_positions = [

            from_si_length(
                x,
                unit_system
            )

            for x in x_positions_si

        ]

        # ----------------------------------------------------
        # Return data
        # ----------------------------------------------------

        return jsonify({

            "success":
                True,

            "vessel_name":
                vessel_name,

            "unit_system":
                unit_system,

            "lpp":
                from_si_length(
                    lpp,
                    unit_system
                ),

            "density":
                from_si_density(
                    density,
                    unit_system
                ),

            "stations":
                number_of_stations,

            "waterlines":
                len(drafts),

            "drafts":
                [
                    from_si_length(
                        x,
                        unit_system
                    )
                    for x in drafts
                ],

            "offsets":
                [
                    [
                        from_si_length(
                            x,
                            unit_system
                        )
                        for x in row
                    ]
                    for row in offsets
                ],

            "x_positions":
                x_positions,

            "results":
                results

        })

    except Exception as e:

        return jsonify({

            "success":
                False,

            "error":
                str(e)

        }), 400


# ============================================================
# RUN LOCAL SERVER
# ============================================================

if __name__ == "__main__":

    app.run(
        host="0.0.0.0",
        port=5000,
        debug=True
    )