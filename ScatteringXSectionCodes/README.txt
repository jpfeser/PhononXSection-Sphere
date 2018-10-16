"TruellXSection_JPFCorrections.m" is for longitudinal phonons, based on Ying and Truell's calculation, but as the name suggests it contains some corrections from the original manuscript (Eq. 20a in particular contains some typos in that paper).
"GetSigma_SphereTrans.m" is for transverse phonons based on calculation is based on: Iwashimazu, Journal of Sound and Vibration, 40(2), 267-271, 1975

The origal algorithm by Ying and Truell (numerical method described by Johnson) is also given for reference.  The results of the two codes (wCorrections and Without) are minor as far as I can tell.

"GetSigmaSphere.m" is a higher level function which requests the scattering cross section (could be either for longitudinal or transverse phonons).
