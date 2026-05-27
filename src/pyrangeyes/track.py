from dataclasses import dataclass, field


@dataclass
class Track:
    """A PyRanges object plus per-track plotting options.

    The wrapped data object is linked, not copied. Options passed here override
    the defaults supplied to :func:`pyrangeyes.plot` for this track.
    """

    data: object
    adapter: str | None = None
    options: dict = field(default_factory=dict)

    def __init__(self, data, adapter=None, **options):
        self.data = data
        self.adapter = adapter
        self.options = dict(options)

    def get(self, key, default=None):
        if key == "adapter":
            return self.adapter if self.adapter is not None else default
        return self.options.get(key, default)

    @property
    def name(self):
        return self.options.get("name")

    def plot(self, **kwargs):
        """Plot this track with :func:`pyrangeyes.plot`.

        Options stored on the track are applied as track-specific options;
        keyword arguments passed here are forwarded as regular plot options.
        """
        from .plot_main import plot

        return plot(self, **kwargs)
