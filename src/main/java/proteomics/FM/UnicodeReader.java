package proteomics.FM;


//package edu.northwestern.at.utils;

/*  Please see the license information at the end of this file. */

import java.io.*;

/**  Generic unicode text reader.
 *
 *  <p>
 *  Uses BOM mark to identify the encoding to be used.
 *  When BOM is not found, uses given default or system encoding.
 *  </p>
 *
 *  <p>
 *  Original pseudocode   : Thomas Weidenfeller<br/>
 *  Implementation tweaked: Aki Nieminen
 *  </p>
 *
 *  <p>
 *  http://www.unicode.org/unicode/faq/utf_bom.html<br/>
 *  BOMs:<br/>
 *    00 00 FE FF    = UTF-32, big-endian<br/>
 *    FF FE 00 00    = UTF-32, little-endian<br/>
 *    EF BB BF       = UTF-8,<br/>
 *    FE FF          = UTF-16, big-endian<br/>
 *    FF FE          = UTF-16, little-endian<br/>
 *  </p>
 */

public class UnicodeReader extends Reader {
    PushbackInputStream internalIn;
    InputStreamReader   internalIn2 = null;
    String              defaultEnc;

    private static final int BOM_SIZE = 4;

    /**  Create Unicode reader.
     *
     * @param in  inputstream to be read
     */

    public UnicodeReader(InputStream in) {
        this( in , null );
    }

    /**  Create Unicode reader.
     *
     * @param in  inputstream to be read
     * @param defaultEnc default encoding if stream does not have
     *                   BOM marker. Give NULL to use system-level default.
     */

    public UnicodeReader(InputStream in, String defaultEnc) {
        internalIn = new PushbackInputStream(in, BOM_SIZE);
        this.defaultEnc = defaultEnc;
    }

    public String getDefaultEncoding() {
        return defaultEnc;
    }

    /**
     * Get stream encoding or NULL if stream is uninitialized.
     * Call init() or read() method to initialize it.
     */
    public String getEncoding() {
        if (internalIn2 == null) return null;
        return internalIn2.getEncoding();
    }

    /**
     * Read-ahead four bytes and check for BOM marks. Extra bytes are
     * unread back to the stream, only BOM bytes are skipped.
     */
    protected void init() throws IOException {
        if (internalIn2 != null) return;

        String encoding;
        byte bom[] = new byte[BOM_SIZE];
        int n, unread;
        n = internalIn.read(bom, 0, bom.length);

        if ( (bom[0] == (byte)0x00) && (bom[1] == (byte)0x00) &&
                (bom[2] == (byte)0xFE) && (bom[3] == (byte)0xFF) ) {
            encoding = "UTF-32BE";
            unread = n - 4;
        } else if ( (bom[0] == (byte)0xFF) && (bom[1] == (byte)0xFE) &&
                (bom[2] == (byte)0x00) && (bom[3] == (byte)0x00) ) {
            encoding = "UTF-32LE";
            unread = n - 4;
        } else if (  (bom[0] == (byte)0xEF) && (bom[1] == (byte)0xBB) &&
                (bom[2] == (byte)0xBF) ) {
            encoding = "UTF-8";
            unread = n - 3;
        } else if ( (bom[0] == (byte)0xFE) && (bom[1] == (byte)0xFF) ) {
            encoding = "UTF-16BE";
            unread = n - 2;
        } else if ( (bom[0] == (byte)0xFF) && (bom[1] == (byte)0xFE) ) {
            encoding = "UTF-16LE";
            unread = n - 2;
        } else {
            // Unicode BOM mark not found, unread all bytes
            encoding = defaultEnc;
            unread = n;
        }
        //System.out.println("read=" + n + ", unread=" + unread);

        if (unread > 0) internalIn.unread(bom, (n - unread), unread);

        // Use given encoding
        if (encoding == null) {
            internalIn2 = new InputStreamReader(internalIn);
        } else {
            internalIn2 = new InputStreamReader(internalIn, encoding);
        }
    }

    public void close() throws IOException {
        init();
        internalIn2.close();
    }

    public int read(char[] cbuf, int off, int len) throws IOException {
        init();
        return internalIn2.read(cbuf, off, len);
    }
}

/*
Copyright (c) 2008, 2009 by Northwestern University.
All rights reserved.

Developed by:
   Academic and Research Technologies
   Northwestern University
   http://www.it.northwestern.edu/about/departments/at/

Permission is hereby granted, free of charge, to any person
obtaining a copy of this software and associated documentation
files (the "Software"), to deal with the Software without
restriction, including without limitation the rights to use,
copy, modify, merge, publish, distribute, sublicense, and/or
sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following
conditions:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimers.

    * Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimers in the documentation and/or other materials provided
      with the distribution.

    * Neither the names of Academic and Research Technologies,
      Northwestern University, nor the names of its contributors may be
      used to endorse or promote products derived from this Software
      without specific prior written permission.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE CONTRIBUTORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS WITH THE SOFTWARE.
*/

