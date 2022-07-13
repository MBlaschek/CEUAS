r"""CostL2 (least squared deviation)"""
from ruptures.costs import NotEnoughPoints

from ruptures.base import BaseCost

import numpy as np


class CostL2(BaseCost):

    r"""
    Least squared deviation.
    """

    model = "l2"

    def __init__(self):
        """Initialize the object."""
        self.signal = None
        self.min_size = 2

    def fit(self, signal) -> "CostL2":
        """Set parameters of the instance.

        Args:
            signal (array): array of shape (n_samples,) or (n_samples, n_features)

        Returns:
            self
        """
        if signal.ndim == 1:
            self.signal = signal.reshape(-1, 1)
        else:
            self.signal = signal

        self.cumsumsquares=np.nancumsum(signal*signal,axis=0)
        self.cumsum=np.nancumsum(signal,axis=0)
        return self

    def error(self, start, end) -> float:
        """Return the approximation cost on the segment [start:end].

        Args:
            start (int): start of the segment
            end (int): end of the segment

        Returns:
            segment cost

        Raises:
            NotEnoughPoints: when the segment is too short (less than `min_size` samples).
        """
        if end - start < self.min_size:
            raise NotEnoughPoints

        #cost=np.nansum(np.nanvar(self.signal[start:end],axis=0) * np.sum(~np.isnan(self.signal[start:end]),axis=0))
        cost2=np.nansum(self.cumsumsquares[end-1,:]-self.cumsumsquares[start,:]-(self.cumsum[end-1,:]-self.cumsum[start,:])**2/(end - start))
        return cost2
