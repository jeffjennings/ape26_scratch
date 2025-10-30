"""
Implementation of APE 26: Removing data storage (representations) from coordinate frames.

This is a prototype that separates reference frames from coordinate data, defining new
classes and simplifying existing ones relative to the existing astropy.coordinates framework.
"""

from __future__ import annotations

from abc import ABC
from dataclasses import dataclass, field, replace
from typing import TYPE_CHECKING, Any, Literal, Union, Self

import numpy as np

from astropy import units as u
from astropy.coordinates import representation as r

from astropy.coordinates.angles import Angle
from astropy.coordinates.baseframe import BaseCoordinateFrame
from astropy.time import Time
from astropy.utils.decorators import lazyproperty
from astropy.coordinates.transformations.graph import TransformGraph
from astropy.coordinates import ICRS, FK5, Galactic

if TYPE_CHECKING:
    from astropy.coordinates.representation import BaseRepresentation, BaseDifferential


@dataclass(frozen=True)
class BaseFrame(ABC):
    """
    Base class for all reference frames without data storage.
    """
    
    def __post_init__(self):
        """Initialize frame attributes after dataclass construction."""
        pass
    
    def __repr__(self):
        """
        String representation of the frame.        
        Returns the frame class name and any frame attributes.
        """
        frame_name = self.__class__.__name__
        
        # Get frame attributes
        if hasattr(self, '__dataclass_fields__'):
            attrs = []
            for field_name in self.__dataclass_fields__:
                value = getattr(self, field_name)
                
                # Format
                if isinstance(value, Time):
                    attrs.append(f"{field_name}={value.value}")
                elif value is not None:
                    attrs.append(f"{field_name}={value}")
            
            if attrs:
                attrs_str = ", ".join(attrs)
                return f"<{frame_name} ({attrs_str})>"
        
        # No attributes to show
        return f"<{frame_name}>"
        
    def transform_data_to(
        self, 
        frame: BaseFrame, 
        data: BaseRepresentation,
        **kwargs
    ) -> BaseRepresentation:
        """
        Transform coordinate data from one frame to another.
        
        Uses the existing astropy transformation graph purely for demonstration purposes, 
        allowing APE 26 frames to leverage all existing transformations.
        
        Parameters
        ----------
        frame : BaseFrame
            The target reference frame to transform to.
        data : BaseRepresentation
            The coordinate data to transform.
            
        Returns
        -------
        BaseRepresentation
            The transformed coordinate data.
        """

        # Check if frames are equivalent (no transformation needed)
        if self.is_equivalent_frame(frame):
            return data
        
        # Convert APE26 frames to legacy frame classes to interface with the current 
        # frame transform graph:
        
        # 1) Get the corresponding legacy frame class of an APE26 frame class 
        # (assumes the former is in global namespace) and add data to it
        fromsys = self._convert_to_legacy_frame(data=data)
        # 2) Get the target legacy frame class
        tosys = self._convert_to_legacy_frame(frame=frame)
        # 3) Use the transformation graph
        try:
            transformed = fromsys.transform_to(tosys)
            return transformed.data
        except Exception as e:
            raise ValueError(
                f"Could not transform from {self.__class__.__name__} to "
                f"{frame.__class__.__name__}: {e}"
            )
    
    def _convert_to_legacy_frame(self, frame=None, data=None) -> BaseCoordinateFrame:
        """
        Convert an APE26 frame to a legacy BaseCoordinateFrame, solely for the purpose
        of using the existing astropy frame transform graph. The 'fromsys' frame will
        be converted to legacy and packed with data; the 'tosys' frame will be converted
        to legacy without data.

        Parameters
        ----------
        frame : BaseFrame, optional
            If provided, treat this as a 'tosys' frame.
        data : BaseRepresentation, optional
            If a 'fromsys' frame, this must be provided.
        
        Returns
        -------
        BaseCoordinateFrame
            The corresponding legacy frame instance.
        """
        if frame is not None:
            return globals()[frame.__class__.__name__.rstrip('Frame')]()
        else:
            assert data is not None, "Data must be provided for 'fromsys' frame conversion."
            legacy_class = globals()[self.__class__.__name__.rstrip('Frame')]
            frame_attrs = self._get_frame_attributes()
            return legacy_class(data, **frame_attrs)

    def _get_frame_attributes(self) -> dict:
        """
        Get frame attributes as a dictionary for legacy frame creation, solely for the 
        purpose of using the existing astropy frame transform graph.
        
        Returns
        -------
        dict
            Dictionary of frame attribute names to values.
        """
        attrs = {}
        if hasattr(self, '__dataclass_fields__'):
            for field_name in self.__dataclass_fields__:
                attrs[field_name] = getattr(self, field_name)
        return attrs
    
    def is_equivalent_frame(self, other: BaseFrame) -> bool:
        """
        Check if this frame is equivalent to another frame.
        
        Parameters
        ----------
        other : BaseFrame
            The other frame to compare with.
            
        Returns
        -------
        bool
            True if frames are equivalent, False otherwise.
        """
        if self.__class__ != other.__class__:
            return False
        
        for field_name in self.__dataclass_fields__:
            self_val = getattr(self, field_name)
            other_val = getattr(other, field_name)
            
            if not self._frameattr_equiv(self_val, other_val):
                return False
        return True
    
    @staticmethod
    def _frameattr_equiv(left_attr, right_attr) -> bool:
        """Check if two frame attributes are equivalent."""
        if left_attr is right_attr:
            return True
        if left_attr is None or right_attr is None:
            return False
        return np.all(left_attr == right_attr)


@dataclass(frozen=True)
class ICRSFrame(BaseFrame):
    """
    International Celestial Reference System (ICRS) frame.
    
    This frame has no specific attributes - it is defined purely by
    the reference frame itself.
    """
    pass

@dataclass(frozen=True)
class GalacticFrame(BaseFrame):
    """
    Galactic coordinate frame.
    
    This frame has no specific attributes.
    """
    pass

@dataclass(frozen=True)
class FK5Frame(BaseFrame):
    """
    FK5 (Fifth Fundamental Catalog) reference frame.
    
    This frame includes an equinox attribute that specifies the epoch
    of the mean equator and equinox.
    """    
    equinox: Time = field(default_factory=lambda: Time("J2000"))

    def __repr__(self):
        """
        String representation matching astropy FK5 format.
        
        Returns
        -------
        str
            Representation like '<FK5Frame (equinox=J2000.000)>'
        """
        # Format equinox in the astropy way
        if self.equinox.scale == 'tt':
            # Use the jyear format for J2000-style epochs
            equinox_str = f"J{self.equinox.jyear:.3f}"
        else:
            equinox_str = str(self.equinox)
        
        return f"<FK5Frame (equinox={equinox_str})>"
    
# =============================================================================
# Coordinate Classes (data + frame)
# =============================================================================

@dataclass
class BaseCoordinate(ABC):
    """
    Base class for data in a reference frame.
    
    This class brings together:
    - A reference frame (BaseFrame instance)
    - Coordinate data (BaseRepresentation instance)
    """
    
    frame: BaseFrame
    data: BaseRepresentation
    
    def __post_init__(self):
        """Validate inputs after dataclass construction."""
        if not isinstance(self.frame, BaseFrame):
            raise TypeError(
                f"frame must be a BaseFrame instance, got {type(self.frame)}"
            )
        if not isinstance(self.data, r.BaseRepresentation):
            raise TypeError(
                f"data must be a BaseRepresentation instance, got {type(self.data)}"
            )
    
    def transform_to(self, new_frame: BaseFrame) -> BaseCoordinate:
        """
        Transform this coordinate to a new reference frame.
        
        Parameters
        ----------
        new_frame : BaseFrame
            The target reference frame.
            
        Returns
        -------
        BaseCoordinate
            A new coordinate instance in the target frame.
        """
        new_data = self.frame.transform_data_to(new_frame, self.data)
        return self.__class__(frame=new_frame, data=new_data)
    
    def separation(self, other: BaseCoordinate) -> Angle:
        """
        Compute on-sky separation to another coordinate.
        
        Parameters
        ----------
        other : BaseCoordinate
            The other coordinate.
            
        Returns
        -------
        Angle
            The angular separation.
        """
        # Transform other to this frame
        other_in_frame = other.transform_to(self.frame)
        
        # Convert to unit spherical for separation calculation
        self_sph = self.data.represent_as(r.UnitSphericalRepresentation)
        other_sph = other_in_frame.data.represent_as(r.UnitSphericalRepresentation)
        
        from astropy.coordinates.angles import angular_separation
        return Angle(
            angular_separation(
                self_sph.lon, self_sph.lat,
                other_sph.lon, other_sph.lat
            ),
            unit=u.degree
        )
    
    def represent_as(
        self, 
        representation_cls: Union[type[BaseRepresentation], str],
        differential_cls: type[BaseDifferential] = None
    ) -> BaseRepresentation:
        """
        Represent the coordinate data in a different representation.
        
        Parameters
        ----------
        representation_cls : BaseRepresentation subclass or str
            The representation class to use, or a string name like 
            'spherical', 'cartesian', 'cylindrical', 'unitspherical'.
        differential_cls : BaseDifferential subclass, optional
            The differential class to use for velocities.
            
        Returns
        -------
        BaseRepresentation
            The data in the new representation.
        """
        # Handle string representation names
        if isinstance(representation_cls, str):
            representation_name = representation_cls.lower()
            
            if representation_name in ('spherical', 'sphericalrepresentation'):
                representation_cls = r.SphericalRepresentation
            elif representation_name in ('cartesian', 'cartesianrepresentation'):
                representation_cls = r.CartesianRepresentation
            elif representation_name in ('cylindrical', 'cylindricalrepresentation'):
                representation_cls = r.CylindricalRepresentation
            elif representation_name in ('unitspherical', 'unitsphericalrepresentation'):
                representation_cls = r.UnitSphericalRepresentation
            elif representation_name in ('physicsspherical', 'physicssphericalrepresentation'):
                representation_cls = r.PhysicsSphericalRepresentation
            else:
                raise ValueError(
                    f"Unknown representation name: {representation_cls}. "
                    "Use 'spherical', 'cartesian', 'cylindrical', or 'unitspherical'."
                )
        
        if differential_cls is not None and 's' in self.data.differentials:
            return self.data.represent_as(representation_cls, differential_cls)
        return self.data.represent_as(representation_cls)
    
    @property
    def shape(self):
        """Shape of the coordinate data."""
        return self.data.shape
    
    @property
    def size(self):
        """Size of the coordinate data."""
        return self.data.size
    
    def __len__(self):
        """Length of the coordinate array."""
        return len(self.data)
    
    def __getitem__(self, item):
        """Index into the coordinate data."""
        return self.__class__(frame=self.frame, data=self.data[item])


@dataclass
class Coordinate(BaseCoordinate):
    """
    Lightweight class for data in a reference frame.
    """
    data: BaseRepresentation
    frame: BaseFrame

    # TODO: remove when astropy supports python 3.13+ (just use copy.replace(self, ...))
    def __replace__(self, **changes) -> Self:
        return replace(self, **changes)
   
    def transform_to(self, new_frame: BaseFrame) -> Self:
        """Transform to a new frame."""
        new_data = self.frame.transform_data_to(new_frame, self.data)
        # TODO: python 3.13+ just use copy.replace(self, ...)
        return replace(self, data=new_data, frame=new_frame)
    
    def __repr__(self):
        """String representation."""
        frame_name = self.frame.__class__.__name__
        data_repr = repr(self.data)
        return f"<Coordinate ({frame_name}): {data_repr}>"


@dataclass
class SkyCoord(BaseCoordinate):
    """
    Feature-rich coordinate class with batteries included.
    
    This class includes:
    - Flexible input parsing
    - Caching of transformations
    - Storage of frame attributes not in current frame
    - Rich API with many convenience methods
    
    This is the primary user-facing coordinate class.
    
    NOTE: This uses a regular (non-frozen) dataclass to support caching.
    """
    
    _cache: dict[str, Any] = field(default_factory=dict, init=False, repr=False)
    _extra_frameattrs: dict[str, Any] = field(
        default_factory=dict, 
        init=False, 
        repr=False
    )
    
    def __init__(
        self, 
        *args,
        frame: Union[BaseFrame, str] = None,
        data: BaseRepresentation = None,
        copy: bool = True,
        **kwargs
    ):
        """
        Initialize SkyCoord with flexible input parsing.
        
        Parameters
        ----------
        *args : positional arguments
            Coordinate components (e.g., ra, dec) or representation.
            Can also be a legacy BaseCoordinateFrame instance.
        frame : BaseFrame or string
            The reference frame. Can be a BaseFrame instance or frame name.
        data : BaseRepresentation
            The coordinate data.
        copy : bool
            Whether to copy the data.
        **kwargs
            Additional frame attributes or coordinate components.
        """
        # Handle single argument that might be a legacy frame or coordinate
        if len(args) == 1 and frame is None and data is None:
            arg = args[0]
            
            # Check if it's a legacy BaseCoordinateFrame
            if isinstance(arg, BaseCoordinateFrame):
                frame = arg.frame
                data = arg.data
                args = ()
            # Check if it's a BaseCoordinate
            elif isinstance(arg, BaseCoordinate):
                frame = arg.frame
                data = arg.data
                args = ()
        
        # Handle other positional arguments (but only if data wasn't already set)
        if args and data is None:
            if len(args) == 1 and isinstance(args[0], r.BaseRepresentation):
                data = args[0]
            elif len(args) >= 2:
                # Assume positional ra, dec (and optionally distance)
                kwargs['ra'] = args[0]
                kwargs['dec'] = args[1]
                if len(args) > 2:
                    kwargs['distance'] = args[2]
        
        # Handle string frame names or None (default to ICRS)
        if isinstance(frame, str) or frame is None:
            if frame is None:
                frame = 'icrs'
            frame = self._frame_from_string(frame, **kwargs)
        
        # Flexible input parsing for coordinate components (only if data is still None)
        if data is None and kwargs:
            data = self._parse_coordinate_data(**kwargs)
        
        # Validate that we have data at this point
        if data is None:
            raise ValueError(
                "Could not determine coordinate data. Please provide coordinate "
                "components (e.g., ra, dec) or a BaseRepresentation object."
            )
        
        if copy and data is not None:
            data = data.copy()
        
        # Initialize using object.__setattr__ for frozen-like behavior
        object.__setattr__(self, 'frame', frame)
        object.__setattr__(self, 'data', data)
        object.__setattr__(self, '_cache', {})
        object.__setattr__(self, '_extra_frameattrs', {})
        
        # Validate using the parent's __post_init__
        BaseCoordinate.__post_init__(self)
    
    def _frame_from_string(self, frame_name: str, **kwargs) -> BaseFrame:
        """Create a frame from a string name."""
        frame_name_lower = frame_name.lower()
        
        if frame_name_lower == "icrs":
            return ICRSFrame()
        elif frame_name_lower == "fk5":
            equinox = kwargs.pop('equinox', None)
            if equinox is not None:
                return FK5Frame(equinox=equinox)
            return FK5Frame()
        elif frame_name_lower == "galactic":
            return GalacticFrame()
        else:
            raise ValueError(f"Unknown frame name: {frame_name}")
    
    def _parse_coordinate_data(self, **kwargs) -> BaseRepresentation:
        """Parse coordinate components into a representation."""
        # Remove representation_type hint if present
        representation_type = kwargs.pop('representation_type', None)
        
        # Parse Cartesian coordinates
        if 'x' in kwargs and 'y' in kwargs and 'z' in kwargs:
            x = kwargs.pop('x')
            y = kwargs.pop('y')
            z = kwargs.pop('z')
            return r.CartesianRepresentation(x=x, y=y, z=z)
        
        # Parse spherical coordinates (ra, dec)
        if 'ra' in kwargs and 'dec' in kwargs:
            ra = kwargs.pop('ra')
            dec = kwargs.pop('dec')
            distance = kwargs.pop('distance', None)
            
            if distance is not None:
                return r.SphericalRepresentation(
                    lon=ra, lat=dec, distance=distance
                )
            else:
                return r.UnitSphericalRepresentation(lon=ra, lat=dec)
        
        # Parse generic lon, lat coordinates
        if 'lon' in kwargs and 'lat' in kwargs:
            lon = kwargs.pop('lon')
            lat = kwargs.pop('lat')
            distance = kwargs.pop('distance', None)
            
            if distance is not None:
                return r.SphericalRepresentation(
                    lon=lon, lat=lat, distance=distance
                )
            else:
                return r.UnitSphericalRepresentation(lon=lon, lat=lat)
        
        # Parse cylindrical coordinates
        if 'rho' in kwargs and 'phi' in kwargs and 'z' in kwargs:
            rho = kwargs.pop('rho')
            phi = kwargs.pop('phi')
            z = kwargs.pop('z')
            return r.CylindricalRepresentation(rho=rho, phi=phi, z=z)
        
        raise ValueError(
            f"Could not parse coordinate data from kwargs: {list(kwargs.keys())}. "
            "Supported formats: (x, y, z), (ra, dec[, distance]), "
            "(lon, lat[, distance]), (rho, phi, z)"
        )
    
    @lazyproperty
    def cartesian(self):
        """Cartesian representation of the coordinates."""
        return self.represent_as(r.CartesianRepresentation)
    
    @lazyproperty
    def spherical(self):
        """Spherical representation of the coordinates."""
        return self.represent_as(r.SphericalRepresentation)

    @lazyproperty
    def cylindrical(self):
        """Cylindrical representation of the coordinates."""
        return self.represent_as(r.CylindricalRepresentation)
    
    def transform_to(self, new_frame: Union[BaseFrame, str]) -> SkyCoord:
        """
        Transform to a new frame, with caching.
        
        Parameters
        ----------
        new_frame : BaseFrame or str
            The target frame.
            
        Returns
        -------
        SkyCoord
            Transformed coordinate.
        """
        if isinstance(new_frame, str):
            new_frame = self._frame_from_string(new_frame)
        
        # Check cache
        cache_key = f'transform_{new_frame.__class__.__name__}_{id(new_frame)}'
        if cache_key in self._cache:
            return self._cache[cache_key]
        
        # Perform transformation using the frame's transform_data_to method
        new_data = self.frame.transform_data_to(new_frame, self.data)
        
        # Create new SkyCoord directly without going through flexible __init__
        result = object.__new__(SkyCoord)
        object.__setattr__(result, 'frame', new_frame)
        object.__setattr__(result, 'data', new_data)
        object.__setattr__(result, '_cache', {})
        object.__setattr__(result, '_extra_frameattrs', {})
        
        # Cache result
        self._cache[cache_key] = result
        return result
    
    def __repr__(self):
        """String representation."""
        frame_name = self.frame.__class__.__name__
        data_repr = repr(self.data)
        return f"<SkyCoord ({frame_name}): {data_repr}>"

